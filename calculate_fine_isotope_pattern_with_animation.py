import pandas as pd
from pyopenms import *
import sys
import IsoSpecPy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from pyopenms import EmpiricalFormula
from scipy.stats import wasserstein_distance
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
import pyopenms as oms
import multiprocessing
import functools

ex_H_dict = {"C":1,"D":1,"E":1,"H":2,"K":2,"N":2,"Q":2,"R":4,"S":1,"T":1,"W":1,"Y":1}

class ObjectiveFunction:
    def __init__(self, formula_dict, observed_masses, observed_intensities, mass_tolerance):
        self.formula_dict = formula_dict
        self.observed_masses = observed_masses
        self.observed_intensities = observed_intensities
        self.mass_tolerance = mass_tolerance
        self.iteration_data = []  # list to store data at each iteration

    def __call__(self, initial_values):
        # Generate theoretical isotope distribution
        isopyobj = build_isopy_obj(initial_values, self.formula_dict)
        isopyobj = isopyobj.binned(bin_size)
        isopyobj.normalize()
        
        theoretical_masses, theoretical_probs = isopyobj.np_masses(), isopyobj.np_probs()
        
        if self.formula_dict['Ion Mode'] == 'neg':
            theoretical_masses -= 1.00727647  # Subtract proton mass
        else:
            theoretical_masses += 1.00727647  # Add proton mass
        # Normalize theoretical probabilities
        theoretical_probs = theoretical_probs / np.max(theoretical_probs)

        # Extract and normalize observed peaks
        extracted_masses, normalized_observed_intensities = extract_peptide_peaks(
            self.observed_masses, self.observed_intensities, theoretical_masses, self.mass_tolerance
        )
        # There was no peptide found in the spectra, return a high number
        if np.sum(normalized_observed_intensities) == 0:
            distance = 100
        else:
            # Calculate the Wasserstein distance
            distance = wasserstein_distance(theoretical_masses, extracted_masses, theoretical_probs, normalized_observed_intensities)

        # Store the data for animation
        self.iteration_data.append({
            'initial_values': initial_values.copy(),
            'theoretical_masses': theoretical_masses.copy(),
            'theoretical_probs': theoretical_probs.copy(),
            'distance': distance,
            'extracted_masses': extracted_masses.copy(),
            'normalized_observed_intensities': normalized_observed_intensities.copy()
        })

        return distance

# Define the function for theoretical digestion
def theoretical_digest(sequence: str, enzyme: str = "Trypsin", missed_cleavages: int = 1, ion_mode: str = "pos"):
    # Initialize the ProteaseDigestion object
    digestion = ProteaseDigestion()
    digestion.setEnzyme(enzyme)
    digestion.setMissedCleavages(missed_cleavages)

    # Perform the digestion
    peptides = []
    digestion.digest(AASequence.fromString(sequence), peptides)

    # Extract the molecular formula for each peptide
    peptide_formulas = []
    for peptide in peptides:
        if peptide.getAverageWeight() < 500 or peptide.getAverageWeight() > 3200:
            continue
        formula = peptide.getFormula()
        peptide_formulas.append({
            "Peptide Sequence": peptide.toString(),
            "Molecular Formula": formula.toString()
        })

    # Convert to DataFrame
    df = pd.DataFrame(peptide_formulas)
    return df

def GaussFilterSpectrum(exp):
    # Modifies exp in place
    gf = GaussFilter()
    gf_params = gf.getParameters()
    gf_params.setValue(b'ppm_tolerance', float(ppm))
    gf_params.setValue(b'use_ppm_tolerance', 'true')
    gf_params.setValue(b'gaussian_width', bin_size)
    gf.setParameters(gf_params)
    gf.filterExperiment(exp)


def FlattenBaseline(exp):
    flatter = MorphologicalFilter()
    flatter.filter(exp)

def GetPeakList(file, signal_to_noise):
    # Create experiment object
    if exp_file.endswith("mzML"):
        exp = MSExperiment()
        MzMLFile().load(file, exp)
    elif exp_file.endswith("mzXML"):
        exp = MSExperiment()
        MzXMLFile().load(file, exp)
    else:
        raise Exception("Error loading file. Must be mzML or mzXML format!")
    GaussFilterSpectrum(exp)
    # Get spectra
    spectra = exp.getSpectra()[0]

    # Flatten the baseline
    FlattenBaseline(spectra)

    # Perform peak picking to get centroid spectra
    centroid_spectra = PeakPicking(spectra, signal_to_noise)
    exp_dummy = MSExperiment()
    exp_dummy.addSpectrum(centroid_spectra)
    centroid_file_name = "centroid.mzML"
    oms.MzMLFile().store(centroid_file_name, exp_dummy)
    return centroid_spectra

def PeakPicking(spectrum, signal_to_noise):
    picker = PeakPickerIterative()
    params = picker.getParameters()
    # Set the signal-to-noise threshold
    params.setValue('signal_to_noise_', signal_to_noise)
    picker.setParameters(params)

    centroid_spec = MSSpectrum()
    picker.pick(spectrum, centroid_spec)

    return centroid_spec

def determine_ex_H(theoretical_df, ion_mode='pos'):


    sequences = theoretical_df['Peptide Sequence']

    peptide_ex_Hs = []
    for sequence in sequences:
        ex_Hs = 2 + len(sequence)
        for letter in sequence:
            if letter in ex_H_dict.keys():
                ex_Hs += ex_H_dict[letter]
        peptide_ex_Hs.append(ex_Hs)

    theoretical_df['Ex. Hydrogens'] = peptide_ex_Hs

    return theoretical_df

def subtract_heavy_atoms(theoretical_df):
    # For each peptide, compute the molecular formula after subtracting heavy atoms

    light_formulas = []
    heavy_C_list = []
    heavy_N_list = []
    heavy_H_list = []

    for idx, row in theoretical_df.iterrows():
        formula_str = row['Molecular Formula']
        ex_H = row['Ex. Hydrogens']

        # Convert formula string to EmpiricalFormula
        formula = EmpiricalFormula(formula_str)

        # Get counts of elements
        element_counts = {}
        for elem in formula.getElementalComposition():
            symbol = elem.decode()
            count = formula.getElementalComposition()[elem]
            element_counts[symbol] = count

        # Compute heavy-labeled hydrogens (non-exchangeable hydrogens)
        H_total = element_counts.get('H', 0)
        H_heavy = H_total - ex_H

        # Number of heavy-labeled carbons and nitrogens (all C and N atoms are heavy-labeled)
        C_heavy = element_counts.get('C', 0)
        N_heavy = element_counts.get('N', 0)

        # Save heavy atom counts
        heavy_C_list.append(C_heavy)
        heavy_N_list.append(N_heavy)
        heavy_H_list.append(H_heavy)

        # Subtract heavy atoms to get light formula
        # H: keep only exchangeable hydrogens
        # +1 for M+H, -1 for M-H
        H_light = ex_H

        # Create new element counts excluding heavy-labeled atoms
        new_element_counts = {}
        for elem, count in element_counts.items():
            if elem == 'C':
                continue  # Exclude all carbons (heavy-labeled)
            elif elem == 'N':
                continue  # Exclude all nitrogens (heavy-labeled)
            elif elem == 'H':
                if H_light > 0:
                    new_element_counts['H'] = H_light
            else:
                new_element_counts[elem] = count  # Keep other elements (O, S, etc.)

        # Create new EmpiricalFormula for the light molecular formula
        new_formula = EmpiricalFormula()
        for elem, count in new_element_counts.items():
            new_formula += EmpiricalFormula(f"{elem}{int(count)}")

        # Add new formula string to list
        light_formulas.append(new_formula.toString())

    # Add new columns to DataFrame
    theoretical_df['Light Molecular Formula'] = light_formulas
    theoretical_df['Heavy C'] = heavy_C_list
    theoretical_df['Heavy H'] = heavy_H_list
    theoretical_df['Heavy N'] = heavy_N_list

    return theoretical_df

def extract_peptide_peaks(observed_masses, observed_intensities, theoretical_masses, mass_tolerance):
    extracted_masses = []
    extracted_intensities = []

    for tm in theoretical_masses:
        # Find observed peaks within the mass tolerance window
        mask = np.abs(observed_masses - tm)/tm * 1e6 <= mass_tolerance
        if np.any(mask):
            # If multiple peaks are within the window, take the one with the highest intensity
            idx = np.where(mask)[0]
            intensities = observed_intensities[idx]
            max_idx = idx[np.argmax(intensities)]
            extracted_masses.append(observed_masses[max_idx])
            extracted_intensities.append(observed_intensities[max_idx])
        else:
            # If no peak is found, append zero intensity
            extracted_masses.append(tm)
            extracted_intensities.append(0.0)

    # Convert to NumPy arrays
    extracted_masses = np.array(extracted_masses)
    extracted_intensities = np.array(extracted_intensities)

    #extracted_masses, extracted_intensities = bin_masses(extracted_masses, extracted_intensities, bin_size)
    #pdb.set_trace()
    # Normalize the intensities so that they sum to 1
    total_intensity = np.sum(extracted_intensities)
    if total_intensity > 0:
        normalized_intensities = extracted_intensities / total_intensity
        normalized_intensities = normalized_intensities / np.max(normalized_intensities)
    else:
        # If total intensity is zero, return the original zero intensities
        normalized_intensities = extracted_intensities

    return extracted_masses, normalized_intensities

def compare_theo_obs_graphically(theo_masses, theo_probs, obs_masses, obs_probs, peptide):

    plt.bar(theo_masses, theo_probs, width=0.03, label='Theoretical Masses')
    plt.bar(obs_masses, obs_probs, width=0.03, label='Observed Masses')
    plt.legend()
    plt.title(f"{peptide} Theo. vs. Exp.")
    plt.xlabel('m/z')
    plt.ylabel("Normalized Intensity")
    plt.savefig(f"{peptide}_Theo_Exp.png")
    plt.close()

    return None

def build_isopy_obj(prob_values, formula_dict):

    prob_h_h = prob_values[0]
    prob_h_n = prob_values[1]
    prob_h_c = prob_values[2]
    prob_l_h = 1 - prob_h_h
    prob_l_n = 1 - prob_h_n
    prob_l_c = 1 - prob_h_c

     # were 0.95 and 0.8 before at seq length 8
    peptide_len = len(formula_dict['Peptide Sequence'])
    if peptide_len < 8:
        prob_coverage = 0.99
    elif peptide_len >= 8 and peptide_len < 12: # Hueristic because the isotope envolope will become long and intersect other envolopes
        prob_coverage = 0.95
    elif peptide_len >= 12 and peptide_len < 16:
        prob_coverage = 0.85
    else:
        prob_coverage = 0.8
    
    isopyobj = IsoSpecPy.IsoTotalProb(prob_coverage, formula_dict['Light Molecular Formula'],
                                        atomCounts=(formula_dict['Heavy C'], formula_dict['Heavy H'], formula_dict['Heavy N']),
                                        isotopeProbabilities=([prob_l_c,prob_h_c],[prob_l_h,prob_h_h],[prob_l_n,prob_h_n]),
                                        isotopeMasses=(isotope_mass["C"],isotope_mass["H"],isotope_mass['N']))

    return isopyobj

def bin_masses(theoretical_masses, theoretical_probs, bin_width):
    min_mass = np.min(theoretical_masses)
    max_mass = np.max(theoretical_masses)
    bins = np.arange(min_mass, max_mass + bin_width, bin_width)
    bin_indices = np.digitize(theoretical_masses, bins)
    binned_masses = []
    binned_probs = []
    for i in range(1, len(bins)):
        indices_in_bin = np.where(bin_indices == i)[0]
        if len(indices_in_bin) > 0:
            masses_in_bin = theoretical_masses[indices_in_bin]
            probs_in_bin = theoretical_probs[indices_in_bin]
            total_prob = np.sum(probs_in_bin)
            if total_prob > 0:
                average_mass = np.sum(masses_in_bin * probs_in_bin) / total_prob
                binned_masses.append(average_mass)
                binned_probs.append(total_prob)
            else:
                # total_prob is zero; skip this bin to avoid NaNs
                continue
    return np.array(binned_masses), np.array(binned_probs)



def fit_isotope_abundances(formula_dict, observed_masses, observed_intensities, initial_guess, ppm):
    # Bounds for isotope abundances (between 0 and 1)
    bounds = [(0.9, 0.999), (0.97,0.999), (0.97,0.999)]
    # Convert initial_guess dict to a list
    initial_values = [initial_guess[element] for element in initial_guess]

    # Compute the initial theoretical masses and probabilities using initial_guess
    isopyobj_init = build_isopy_obj(initial_values, formula_dict)
    isopyobj_init = isopyobj_init.binned(bin_size)
    isopyobj_init.normalize()

    theoretical_masses_init, theoretical_probs_init = isopyobj_init.np_masses(), isopyobj_init.np_probs()
    
    if formula_dict['Ion Mode'] == 'neg':
        theoretical_masses_init -= 1.00727647  # Subtract proton mass
    else:
        theoretical_masses_init += 1.00727647  # Add proton mass

    theoretical_probs_init = theoretical_probs_init / np.max(theoretical_probs_init)

    # Extract observed peaks using initial theoretical masses
    extracted_masses_fixed, normalized_observed_intensities_fixed = extract_peptide_peaks(
        observed_masses, observed_intensities, theoretical_masses_init, ppm
    )

    # Create instance of ObjectiveFunction
    obj_func = ObjectiveFunction(formula_dict, observed_masses, observed_intensities, ppm)

    # Optimization
    result = minimize(
        obj_func,
        initial_values,
        bounds=bounds,
        method='L-BFGS-B',
        options={'maxfun':50}
    )

    # Map the result back to isotope abundances
    optimized_abundances = dict(zip(initial_guess.keys(), result.x))

    # Compute Wasserstein distance between final theoretical and observed intensities
    final_data = obj_func.iteration_data[-1]
    final_theoretical_masses = final_data['theoretical_masses']
    final_theoretical_probs = final_data['theoretical_probs']
    final_extracted_masses = final_data['extracted_masses']
    final_normalized_observed_intensities = final_data['normalized_observed_intensities']

    if np.sum(final_normalized_observed_intensities) == 0:
        wasserstein_dist = 100
    else:
        wasserstein_dist = wasserstein_distance(
            final_theoretical_masses, final_extracted_masses,
            final_theoretical_probs, final_normalized_observed_intensities
        )

    # Now, we can create the animation using 'obj_func.iteration_data' and the fixed observed masses and intensities
    create_animation(
        obj_func.iteration_data,
        formula_dict['Peptide Sequence'],
        extracted_masses_fixed,
        normalized_observed_intensities_fixed
    )

    return optimized_abundances, wasserstein_dist, obj_func.iteration_data, extracted_masses_fixed, normalized_observed_intensities_fixed

def create_animation(iteration_data, peptide_sequence, observed_masses_fixed, normalized_observed_intensities_fixed):
    if len(iteration_data) < 2:
        print("Not enough data to create animation.")
        return

    # Extract the error values
    errors = [data['distance'] for data in iteration_data]
    max_error = max(errors)
    y_max_error = max_error * 1.1  # Add some margin

    num_iterations = len(iteration_data)
    x_iterations = range(num_iterations)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # Get the theoretical masses and probabilities from the last frame
    data = iteration_data[-1]
    theoretical_masses = data['theoretical_masses']
    theoretical_probs = data['theoretical_probs']

    # Find the mass with the maximum theoretical probability
    peak_index = np.argmax(theoretical_probs)
    peak_mz = theoretical_masses[peak_index]

    # Set x-axis limits
    x_min = peak_mz - 8
    x_max = peak_mz + 3

    # Set up the first subplot (spectra)
    ax1.set_title(f"Optimization for peptide {peptide_sequence}")
    ax1.set_xlabel('m/z')
    ax1.set_ylabel('Normalized Intensity')
    ax1.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    ax1.set_xlim(x_min, x_max)

    # Set up the second subplot (error)
    ax2.set_title('Wasserstein Distance')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Error')
    ax2.set_ylim(0, y_max_error)
    ax2.set_xlim(0, num_iterations)
    ax2.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    def init():
        return []

    def update(frame):
        data = iteration_data[frame]
        theoretical_masses = data['theoretical_masses']
        theoretical_probs = data['theoretical_probs']
        current_error = data['distance']

        # Update the spectra subplot
        ax1.clear()


        ax1.bar(observed_masses_fixed, normalized_observed_intensities_fixed, width=0.3, label='Observed Masses',
                color='#8ACE00', edgecolor='black')  # Set color to blue
        markerline, stemlines, baseline = ax1.stem(
            theoretical_masses,
            theoretical_probs,
            linefmt='black',
            markerfmt='black',
            label='Theoretical'
        )

        plt.setp(markerline, 'markersize', 4)
        plt.setp(baseline,alpha=0)
        ax1.set_title(f"Optimization for peptide {peptide_sequence}")
        ax1.set_xlabel('m/z')
        ax1.set_ylabel('Normalized Intensity')
        ax1.legend()
        ax1.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(0, 1)

        # Update the error subplot
        ax2.clear()
        # Plot the error values up to the current frame without markers
        ax2.plot(x_iterations[:frame+1], errors[:frame+1], color='blue')
        ax2.set_title('Wasserstein Distance')
        ax2.set_xlabel('Iteration')
        ax2.set_ylabel('Error')
        ax2.set_ylim(0, y_max_error)
        ax2.set_xlim(0, num_iterations)
        ax2.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

        return []

    # Set the desired frames per second (increase to speed up the animation)
    fps = 8  # Adjust this value as needed

    # Create the animation with a smaller interval between frames
    ani = animation.FuncAnimation(
        fig, update, frames=num_iterations,
        init_func=init, blit=False, interval=50  # Decrease interval to speed up
    )

    # Create a PillowWriter instance with the desired fps
    writer = PillowWriter(fps=fps)

    # Save the animation with the specified writer
    ani.save(f"{peptide_sequence}_optimization_animation.gif", writer=writer)
    plt.close()

def process_peptide(row_dict, mode, obs_peak_values, obs_peak_intensities, initial_guess, ppm):
    peptide_dict = row_dict.copy()
    peptide_dict['Ion Mode'] = mode 
    optimized_abundances, wasserstein_dist, iteration_data, extracted_masses_fixed, normalized_observed_intensities_fixed = fit_isotope_abundances(
        peptide_dict, obs_peak_values, obs_peak_intensities, initial_guess, ppm
    )

    print(f"Optimized Isotope Abundances for Peptide {peptide_dict['Peptide Sequence']}: {optimized_abundances}")
    print(f"Wasserstein Distance for Peptide {peptide_dict['Peptide Sequence']}: {wasserstein_dist}")

    # Return data
    return {
        'peptide_sequence': peptide_dict['Peptide Sequence'],
        'wasserstein_distance': wasserstein_dist,
        'optimized_abundances': optimized_abundances,  # Store optimized abundances here
        'iteration_data': iteration_data,
        'observed_masses_fixed': extracted_masses_fixed,
        'normalized_observed_intensities_fixed': normalized_observed_intensities_fixed
    }


ebd = ElementDB()
hydrogen = ebd.getElement("H")
oxygen = ebd.getElement("O")
sulfur = ebd.getElement("S")
carbon = ebd.getElement("C")
nitrogen = ebd.getElement("N")

# Define natural abundances
isotope_dis = {
    "H": [i.getIntensity() for i in hydrogen.getIsotopeDistribution().getContainer()][0:2],
    "O": [i.getIntensity() for i in oxygen.getIsotopeDistribution().getContainer()],
    "S": [i.getIntensity() for i in sulfur.getIsotopeDistribution().getContainer()],
    "C": [i.getIntensity() for i in carbon.getIsotopeDistribution().getContainer()],
    "N": [i.getIntensity() for i in nitrogen.getIsotopeDistribution().getContainer()],
}

# Define masses of isotopes
isotope_mass = {
    "H": [i.getMZ() for i in hydrogen.getIsotopeDistribution().getContainer()][0:2],
    "O": [i.getMZ() for i in oxygen.getIsotopeDistribution().getContainer()],
    "S": [i.getMZ() for i in sulfur.getIsotopeDistribution().getContainer()],
    "C": [i.getMZ() for i in carbon.getIsotopeDistribution().getContainer()],
    "N": [i.getMZ() for i in nitrogen.getIsotopeDistribution().getContainer()],
}

exp_file = sys.argv[1]
enzyme = sys.argv[2]
if enzyme.lower() == 'gluc':
    enzyme = 'glutamyl endopeptidase'
ppm = int(sys.argv[3])
mode = sys.argv[4] 
bin_size = float(sys.argv[5])

# Example sequence
sequence = "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"

# Perform the theoretical digestion
theoretical_df = theoretical_digest(sequence, enzyme, ion_mode=mode)
theoretical_df = theoretical_df.drop_duplicates(subset=['Peptide Sequence'])

# Find exchangeable protons
theoretical_df = determine_ex_H(theoretical_df, mode)

# Decompose the formulas for IsoSpecPy
theoretical_df = subtract_heavy_atoms(theoretical_df)

initial_guess = {"H": 0.94, "N": 0.989, "C": 0.989}

exp_spectrum = GetPeakList(exp_file, 1.5)
obs_peak_values, obs_peak_intensities = exp_spectrum.get_peaks()

# List to collect data for all peptides
peptides_data = []

# Print peptides to be processed
print("Peptides to be processed:")
print(theoretical_df['Peptide Sequence'].tolist())

peptides_data = []

# Prepare the data
rows = theoretical_df.to_dict('records')

# Use functools.partial to pass additional arguments to the worker function
process_peptide_partial = functools.partial(
    process_peptide,
    mode=mode,
    obs_peak_values=obs_peak_values,
    obs_peak_intensities=obs_peak_intensities,
    initial_guess=initial_guess,
    ppm=ppm
)

with multiprocessing.Pool() as pool:
    results = pool.map(process_peptide_partial, rows)
peptides_data = results

# Sort peptides by Wasserstein distance and select the 4 with the lowest distances
sorted_peptides = sorted(peptides_data, key=lambda x: x['wasserstein_distance'])

csv_data = []
for peptide in sorted_peptides:
    optimized_abundances = peptide['optimized_abundances']
    data = {
        'Peptide Sequence': peptide['peptide_sequence'],
        'Wasserstein Distance': peptide['wasserstein_distance'],
        'Optimized H Abundance': optimized_abundances['H'],
        'Optimized N Abundance': optimized_abundances['N'],
        'Optimized C Abundance': optimized_abundances['C']
    }
    csv_data.append(data)

# Convert to DataFrame
df = pd.DataFrame(csv_data)

# Save to CSV
df.to_csv('sorted_peptides.csv', index=False)
print("Sorted peptides have been saved to 'sorted_peptides.csv'")


###### NMR DATA #########
def find_peptide_positions(peptide, full_sequence):
    positions = []
    peptide = peptide.strip().upper()
    full_sequence = full_sequence.upper()
    index = full_sequence.find(peptide)
    while index != -1:
        # Positions are 1-based indexing
        start_pos = index + 1
        end_pos = start_pos + len(peptide) - 1
        positions.extend(range(start_pos, end_pos + 1))
        # Look for next occurrence
        index = full_sequence.find(peptide, index + 1)
    return positions

# Define the number of side-chain hydrogens beyond the CB carbonn
side_chain_dueteriums = {
    'A': 0,  # Alanine
    'C': 0,  # Cysteine
    'D': 0,  # Aspartic acid
    'E': 2,  # Glutamic acid
    'F': 5,  # Phenylalanine
    'G': 0,  # Glycine
    'H': 2,  # Histidine
    'I': 8,  # Isoleucine
    'K': 6,  # Lysine
    'L': 7,  # Leucine
    'M': 5,  # Methionine
    'N': 0,  # Asparagine
    'P': 4,  # Proline
    'Q': 2,  # Glutamine
    'R': 4,  # Arginine
    'S': 0,  # Serine
    'T': 3,  # Threonine
    'V': 6,  # Valine
    'W': 5,  # Tryptophan
    'Y': 4   # Tyrosine
}

cb_dueteriums = {
    'A': 3,  # Alanine
    'C': 2,  # Cysteine
    'D': 2,  # Aspartic acid
    'E': 2,  # Glutamic acid
    'F': 2,  # Phenylalanine
    'G': 0,  # Glycine
    'H': 2,  # Histidine
    'I': 1,  # Isoleucine
    'K': 2,  # Lysine
    'L': 1,  # Leucine
    'M': 2,  # Methionine
    'N': 2,  # Asparagine
    'P': 2,  # Proline
    'Q': 2,  # Glutamine
    'R': 2,  # Arginine
    'S': 2,  # Serine
    'T': 1,  # Threonine
    'V': 1,  # Valine
    'W': 2,  # Tryptophan
    'Y': 2   # Tyrosine

}

ca_dueteriums = {
    'A': 1,  # Alanine
    'C': 1,  # Cysteine
    'D': 1,  # Aspartic acid
    'E': 1,  # Glutamic acid
    'F': 1,  # Phenylalanine
    'G': 2,  # Glycine
    'H': 1,  # Histidine
    'I': 1,  # Isoleucine
    'K': 1,  # Lysine
    'L': 1,  # Leucine
    'M': 1,  # Methionine
    'N': 1,  # Asparagine
    'P': 1,  # Proline
    'Q': 1,  # Glutamine
    'R': 1,  # Arginine
    'S': 1,  # Serine
    'T': 1,  # Threonine
    'V': 1,  # Valine
    'W': 1,  # Tryptophan
    'Y': 1   # Tyrosine

}

# Step 1: Read the CA and CB dataframes
ca_df = pd.read_csv('Asyn_deuteration_fraction_CA.txt', delimiter=' ')
cb_df = pd.read_csv('Asyn_deuteration_fraction_CB.txt', delimiter=' ')

# Step 2: Impute missing data
ca_df['fracD'].fillna(ca_df['fracD'].mean(), inplace=True)
cb_df['fracD'].fillna(cb_df['fracD'].mean(), inplace=True)

# Extract amino acid code from 'residue' column
ca_df['amino_acid'] = ca_df['residue'].str[0]
cb_df['amino_acid'] = cb_df['residue'].str[0]

ca_df['res_num'] = ca_df['residue'].str[1:].astype(int)
cb_df['res_num'] = cb_df['residue'].str[1:].astype(int)

# Step 3: Calculate average deuteration fractions for CA and CB per amino acid
ca_avg = ca_df.groupby('amino_acid')['fracD'].mean().reset_index()
cb_avg = cb_df.groupby('amino_acid')['fracD'].mean().reset_index()

# Merge the averages
avg_df = pd.merge(ca_avg, cb_avg, on='amino_acid', how='outer', suffixes=('_ca', '_cb'))

# Fill any missing values after merge
avg_df['fracD_ca'].fillna(avg_df['fracD_ca'].mean(), inplace=True)
avg_df['fracD_cb'].fillna(avg_df['fracD_cb'].mean(), inplace=True)

# Compute fracD_avg as the weighted average of fracD_ca and fracD_cb
for _, (aa, fracD_ca, fracD_cb) in avg_df.iterrows():
    avg_df['fracD_avg'] = (avg_df['fracD_ca'] * ca_dueteriums[aa] + avg_df['fracD_cb'] * cb_dueteriums[aa]) / (ca_dueteriums[aa] + cb_dueteriums[aa])

# Overall average deuteration fraction
no_g_no_n = avg_df[~avg_df['amino_acid'].isin(['G', 'N', 'T'])]
overall_fracD_avg = no_g_no_n['fracD_avg'].mean()

# Initialize a dictionary to hold the updated average deuteration fractions
aa_total_fracD = {}
avg_df = avg_df.set_index('amino_acid')

for aa in avg_df.index:
    # Total number of carbons: CA + CB + side-chain carbons
    total_d = ca_dueteriums[aa] + cb_dueteriums[aa] + side_chain_dueteriums[aa]
    # Weighted average by number of protons
    percent_d_residue = (ca_dueteriums[aa] * avg_df.loc[aa]['fracD_ca'] + cb_dueteriums[aa] * avg_df.loc[aa]['fracD_cb'] + side_chain_dueteriums[aa] * overall_fracD_avg)/total_d
    # Add the weighted average for the per residue
    aa_total_fracD[aa] = percent_d_residue

# Step 5: Compute theoretical deuteration percentage for each peptide
theoretical_percentages = []
ca_df = ca_df.set_index('res_num')
cb_df = cb_df.set_index('res_num')

# Initialize lists to store results
theoretical_percentages = []
num_imputations_list = []
total_h_list = []

for idx, peptide in theoretical_df.iterrows():
    peptide_sequence = peptide['Peptide Sequence']
    total_d_per_list = []
    total_d_list = []
    num_imputations = 0  # Initialize counter for this peptide
    total_analyzed_hydrogens = 0 # For calculating the percent imputations

    sequence_positions = find_peptide_positions(peptide_sequence, sequence)

    for aa_idx in sequence_positions:
        aa = sequence[aa_idx - 1]  # Adjust for 0-based indexing

        # Use the experimental values where possible
        if aa_idx in ca_df.index:
            fracD_ca = ca_df.loc[aa_idx]['fracD']
        else:
            fracD_ca = avg_df.loc[aa]['fracD_ca']
            num_imputations += ca_dueteriums[aa]  # Increment counter when imputing 'fracD_ca'
        total_analyzed_hydrogens += ca_dueteriums[aa]

        if aa_idx in cb_df.index:
            fracD_cb = cb_df.loc[aa_idx]['fracD']
        else:
            fracD_cb = avg_df.loc[aa]['fracD_cb']
            num_imputations += cb_dueteriums[aa]  # Increment counter when imputing 'fracD_cb'
        total_analyzed_hydrogens += cb_dueteriums[aa]

        # Number of exchangeable hydrogens for CA, CB, and side chains
        n_ca_d = ca_dueteriums[aa]
        n_cb_d = cb_dueteriums[aa]
        n_side_d = side_chain_dueteriums[aa]

        # Total deuterated hydrogens for this residue
        total_d_per = (fracD_ca * n_ca_d) + (fracD_cb * n_cb_d) + (overall_fracD_avg * n_side_d)

        # Total exchangeable hydrogens for this residue
        total_d = n_ca_d + n_cb_d + n_side_d
        num_imputations += n_side_d
        total_analyzed_hydrogens += n_side_d

        total_d_per_list.append(total_d_per)
        total_d_list.append(total_d)

    # Sum over all residues to get total deuterated hydrogens and total exchangeable hydrogens
    total_deuterated = sum(total_d_per_list)
    total_exchangeable = sum(total_d_list)

    # Compute overall deuteration fraction for the peptide
    avg_fracD_peptide = total_deuterated / total_exchangeable

    # Compute theoretical deuteration percentage
    theoretical_deuteration_percentage = avg_fracD_peptide * 100

    # Append results to lists
    theoretical_percentages.append(theoretical_deuteration_percentage)
    num_imputations_list.append(num_imputations)
    total_h_list.append(total_analyzed_hydrogens)

# Add theoretical deuteration percentages and number of imputations to 'theoretical_df'
theoretical_df['Theoretical_Deuteration_%'] = theoretical_percentages
theoretical_df['Num_Imputations'] = num_imputations_list

# Calculate the total number of possible data points per peptide (each residue has 'fracD_ca' and 'fracD_cb')
total_data_points_list = [2 * len(find_peptide_positions(peptide['Peptide Sequence'], sequence)) for _, peptide in theoretical_df.iterrows()]

# Calculate the percentage of imputed values per peptide
imputation_percentage = np.array(num_imputations_list) / np.array(total_h_list)

# Add imputation percentage to the dataframe
theoretical_df['Imputation_%'] = imputation_percentage * 100

# Add theoretical deuteration percentages to theoretical_df
theoretical_df['Theoretical_Deuteration_%'] = theoretical_percentages

# Display the updated theoretical_df
print(theoretical_df)

csv_data_df = pd.DataFrame(csv_data)

# Step 2: Ensure the peptide sequence columns are aligned for merging
# In theoretical_df, the sequences are under 'Sequence'; we'll rename it to 'Peptide Sequence'
theoretical_df = theoretical_df.rename(columns={'Sequence': 'Peptide Sequence'})

# Step 3: Merge theoretical_df with csv_data_df on 'Peptide Sequence'
merged_df = pd.merge(theoretical_df, csv_data_df, on='Peptide Sequence', how='left')

# Step 4: Compute 'Mass_Spec_Deuteration_%' = (1 - 'Optimized H Abundance') * 100
merged_df['Mass_Spec_Deuteration_%'] = (merged_df['Optimized H Abundance']) * 100

# Step 5: Compute the difference between theoretical and mass spec deuteration percentages
merged_df['Difference'] = merged_df['Theoretical_Deuteration_%'] - merged_df['Mass_Spec_Deuteration_%']

# Filtering out probably not real results
filtered_df = merged_df[(merged_df['Wasserstein Distance'] > 0.05) & (merged_df['Wasserstein Distance'] < 0.4)]
filtered_df = filtered_df[np.abs(filtered_df['Difference']) < 0.5]
# Display the updated theoretical_df
print(filtered_df[['Peptide Sequence', 'Theoretical_Deuteration_%', 'Mass_Spec_Deuteration_%', 'Difference']])

file_name = exp_file.split(".")[0] + "_processed.csv"
filtered_df.to_csv(file_name)

