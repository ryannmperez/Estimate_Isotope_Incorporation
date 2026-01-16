import pandas as pd
from pyopenms import *
import sys
import IsoSpecPy
import numpy as np
import matplotlib.pyplot as plt
from pyopenms import EmpiricalFormula
import pyopenms as oms
import multiprocessing
import functools

ex_H_dict = {"C":1,"D":1,"E":1,"H":2,"K":2,"N":2,"Q":2,"R":4,"S":1,"T":1,"W":1,"Y":1}

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

    # Normalize the intensities so that they sum to 1
    total_intensity = np.sum(extracted_intensities)
    if total_intensity > 0:
        normalized_intensities = extracted_intensities / total_intensity
        normalized_intensities = normalized_intensities / np.max(normalized_intensities)
    else:
        # If total intensity is zero, return the original zero intensities
        normalized_intensities = extracted_intensities

    return extracted_masses, normalized_intensities

def build_isopy_obj(prob_values, formula_dict):

    prob_h_h = prob_values[0]
    prob_h_n = 0.99
    prob_h_c = 0.99
    prob_l_h = 1 - prob_h_h
    prob_l_n = 1 - prob_h_n
    prob_l_c = 1 - prob_h_c

    peptide_len = len(formula_dict['Peptide Sequence'])
    if peptide_len < 8:
        prob_coverage = 0.99
    elif peptide_len >= 8 and peptide_len < 12:
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

def fit_isotope_abundances(formula_dict, observed_masses, observed_intensities, ppm):
    percentages = np.arange(0.94, 0.9801, 0.001)  # From 0.94 to 0.98, step 0.001
    best_distance = float('inf')
    best_percentage = None
    best_theoretical_masses = None
    best_theoretical_probs = None
    best_extracted_masses = None
    best_normalized_observed_intensities = None

    for percentage in percentages:
        prob_values = [percentage]  # Only H percentage

        # Generate theoretical isotope distribution
        isopyobj = build_isopy_obj(prob_values, formula_dict)
        isopyobj = isopyobj.binned(bin_size)
        isopyobj.normalize()
        theoretical_masses, theoretical_probs = isopyobj.np_masses(), isopyobj.np_probs()
        # Adjust masses for ion mode
        if formula_dict['Ion Mode'] == 'neg':
            theoretical_masses -= 1.00727647  # Subtract proton mass
        else:
            theoretical_masses += 1.00727647  # Add proton mass
        # Normalize theoretical probabilities
        theoretical_probs = theoretical_probs / np.max(theoretical_probs)

        # Extract observed peaks
        extracted_masses, normalized_observed_intensities = extract_peptide_peaks(
            observed_masses, observed_intensities, theoretical_masses, ppm
        )
        # Compute distance
        if np.sum(normalized_observed_intensities) == 0:
            distance = 100
        else:
            # Find the index of the theoretical maximum
            theo_max_idx = np.argmax(theoretical_probs)
            theo_max_mz = theoretical_masses[theo_max_idx]
            # Find the index of the observed maximum intensity
            obs_max_idx = np.argmax(normalized_observed_intensities)
            obs_max_mz = extracted_masses[obs_max_idx]
            # Check if theoretical maximum m/z and observed maximum m/z are close (within mass tolerance)
            mass_difference_ppm = np.abs(theo_max_mz - obs_max_mz) / theo_max_mz * 1e6
            if mass_difference_ppm <= ppm:
                # They are close, calculate the difference in intensities at the maximum m/z
                observed_intensity_at_theo_max = normalized_observed_intensities[theo_max_idx]
                theoretical_intensity_at_theo_max = theoretical_probs[theo_max_idx]
                # Compute the absolute difference
                distance = abs(theoretical_intensity_at_theo_max - observed_intensity_at_theo_max)
            else:
                # Not close, set high distance
                distance = 100
        # Update best percentage if distance is better
        if distance < best_distance:
            best_distance = distance
            best_percentage = percentage
            best_theoretical_masses = theoretical_masses
            best_theoretical_probs = theoretical_probs
            best_extracted_masses = extracted_masses
            best_normalized_observed_intensities = normalized_observed_intensities

    # Return best percentage and associated data
    optimized_abundances = {'H': best_percentage}

    # Recompute the theoretical masses and probabilities using the optimized parameters
    isopyobj_opt = build_isopy_obj([best_percentage], formula_dict)
    isopyobj_opt.normalize()
    theoretical_masses_opt, theoretical_probs_opt = isopyobj_opt.np_masses(), isopyobj_opt.np_probs()
    if formula_dict['Ion Mode'] == 'neg':
        theoretical_masses_opt -= 1.00727647  # Subtract proton mass
    else:
        theoretical_masses_opt += 1.00727647  # Add proton mass

    # Normalize theoretical probabilities
    theoretical_probs_opt = theoretical_probs_opt / np.max(theoretical_probs_opt)

    # Extract the observed spectrum in the relevant mass range for plotting
    mass_min = theoretical_masses_opt.min() - 1.0  # Extend range as needed
    mass_max = theoretical_masses_opt.max() + 1.0
    observed_indices = np.where((observed_masses >= mass_min) & (observed_masses <= mass_max))
    observed_masses_in_range = observed_masses[observed_indices]
    observed_intensities_in_range = observed_intensities[observed_indices]
    # Normalize observed intensities
    if observed_intensities_in_range.sum() > 0:
        observed_intensities_in_range = observed_intensities_in_range / observed_intensities_in_range.max()
    else:
        observed_intensities_in_range = observed_intensities_in_range

    return optimized_abundances, best_distance, theoretical_masses_opt, theoretical_probs_opt, observed_masses_in_range, observed_intensities_in_range

def process_peptide(row_dict, mode, obs_peak_values, obs_peak_intensities, ppm):
    peptide_dict = row_dict.copy()
    peptide_dict['Ion Mode'] = mode 
    optimized_abundances, wasserstein_dist, theoretical_masses_opt, theoretical_probs_opt, observed_masses_in_range, observed_intensities_in_range = fit_isotope_abundances(
        peptide_dict, obs_peak_values, obs_peak_intensities, ppm
    )

    print(f"Optimized Isotope Abundances for Peptide {peptide_dict['Peptide Sequence']}: {optimized_abundances}")
    print(f"Distance for Peptide {peptide_dict['Peptide Sequence']}: {wasserstein_dist}")

    # Return data
    return {
        'peptide_sequence': peptide_dict['Peptide Sequence'],
        'distance': wasserstein_dist,
        'optimized_abundances': optimized_abundances,
        'theoretical_masses': theoretical_masses_opt,
        'theoretical_intensities': theoretical_probs_opt,
        'observed_masses': observed_masses_in_range,
        'observed_intensities': observed_intensities_in_range
    }

def plot_spectra(peptide_data):
    peptide_sequence = peptide_data['peptide_sequence']
    theoretical_masses = peptide_data['theoretical_masses']
    theoretical_intensities = peptide_data['theoretical_intensities']
    observed_masses = peptide_data['observed_masses']
    observed_intensities = peptide_data['observed_intensities']

    plt.figure(figsize=(10, 6))
    plt.title(f"Peptide: {peptide_sequence}")
    plt.xlabel('m/z')
    plt.ylabel('Normalized Intensity')

    # Plot observed spectrum
    plt.bar(observed_masses, observed_intensities, width=0.1, color='gray', label='Observed Spectrum', alpha=0.7)

    # Plot theoretical isotope pattern
    plt.stem(theoretical_masses, theoretical_intensities, linefmt='r-', markerfmt='ro', basefmt=' ', label='Theoretical Isotope Pattern')

    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{peptide_sequence}_spectrum.png")
    plt.close()

    print(f"Saved plot for peptide {peptide_sequence} as {peptide_sequence}_spectrum.png")

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

if __name__ == '__main__':
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

    exp_spectrum = GetPeakList(exp_file, 1.5)
    obs_peak_values, obs_peak_intensities = exp_spectrum.get_peaks()

    # List to collect data for all peptides
    peptides_data = []

    # Print peptides to be processed
    print("Peptides to be processed:")
    print(theoretical_df['Peptide Sequence'].tolist())

    # Prepare the data
    rows = theoretical_df.to_dict('records')

    # Use functools.partial to pass additional arguments to the worker function
    process_peptide_partial = functools.partial(
        process_peptide,
        mode=mode,
        obs_peak_values=obs_peak_values,
        obs_peak_intensities=obs_peak_intensities,
        ppm=ppm
    )

    with multiprocessing.Pool() as pool:
        results = pool.map(process_peptide_partial, rows)
    peptides_data = results

    # Sort peptides by distance and select the ones with the lowest distances
    sorted_peptides = sorted(peptides_data, key=lambda x: x['distance'])

    csv_data = []
    for peptide in sorted_peptides:
        optimized_abundances = peptide['optimized_abundances']
        data = {
            'Peptide Sequence': peptide['peptide_sequence'],
            'Distance': peptide['distance'],
            'Optimized H Abundance': optimized_abundances['H']
        }
        csv_data.append(data)
        # Plot the spectra
        plot_spectra(peptide)

    # Convert to DataFrame
    df = pd.DataFrame(csv_data)

    # Save to CSV
    df.to_csv('sorted_peptides.csv', index=False)
    print("Sorted peptides have been saved to 'sorted_peptides.csv'")

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
    filtered_df = merged_df[np.abs(merged_df['Difference']) < 0.5]
    # Display the updated theoretical_df
    print(filtered_df[['Peptide Sequence', 'Theoretical_Deuteration_%', 'Mass_Spec_Deuteration_%', 'Difference']])

    file_name = exp_file.split(".")[0] + "_processed.csv"
    filtered_df.to_csv(file_name)
