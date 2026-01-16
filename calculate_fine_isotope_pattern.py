import pandas as pd
from pyopenms import *
import sys
import IsoSpecPy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from pyopenms import EmpiricalFormula
from scipy.stats import wasserstein_distance
import pdb

# Define the function for theoretical digestion
def theoretical_digest(sequence: str, enzyme: str = "Trypsin", missed_cleavages: int = 1):
    # Initialize the ProteaseDigestion object
    digestion = ProteaseDigestion()
    
    digestion.setEnzyme("glutamyl endopeptidase")
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

    gf = SavitzkyGolayFilter()
    gf_params = gf.getParameters()
    #gf_params.setValue('frame_length', 0.2)
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
        raise "Error loading file. Must be MzML or MzXML format!"
    GaussFilterSpectrum(exp)
    # Get spectra
    spectra = exp.getSpectra()[0]
    
    # Apply Gaussian filter to smooth spectra

    
    # Flatten the baseline
    FlattenBaseline(spectra)
    
    # Perform peak picking to get centroid spectra
    centroid_spectra = PeakPicking(spectra, signal_to_noise)
    
    return centroid_spectra

def PeakPicking(spectrum, signal_to_noise):
    picker = PeakPickerIterative()
    params = picker.getParameters()
    # Set the signal-to-noise threshold
    params.setValue('signal_to_noise', signal_to_noise)
    picker.setParameters(params)
    
    centroid_spec = MSSpectrum()
    picker.pick(spectrum, centroid_spec)
    
    return centroid_spec

def determine_ex_H(theoretical_df):

    ex_H_dict = {"C":1,"D":1,"E":1,"H":2,"K":2,"N":2,"Q":2,"R":4,"S":1,"T":1,"W":1,"Y":1}
    
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

import matplotlib.pyplot as plt

def compare_theo_obs_graphically(theo_masses, theo_probs, obs_masses, obs_probs, peptide):
    plt.figure(figsize=(10, 6))

    # Plot theoretical masses with transparency
    plt.bar(theo_masses, theo_probs, width=0.1, label='Theoretical Masses',
            alpha=1, color='navy', edgecolor='black')

    # Plot observed masses with transparency
    plt.bar(obs_masses, obs_probs, width=0.1, label='Observed Masses',
            alpha=0.5, color='salmon', edgecolor='black')

    # Add gridlines for better readability
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # Customize title and labels
    plt.title(f"{peptide} Theoretical vs. Experimental Masses", fontsize=14)
    plt.xlabel('m/z', fontsize=12)
    plt.ylabel("Normalized Intensity", fontsize=12)

    # Add legend with a fancy box
    plt.legend(frameon=True, fancybox=True, shadow=True)

    # Save the plot as a PNG file
    plt.savefig(f"{peptide}_Theo_Exp.png", dpi=300, bbox_inches='tight')
    plt.close()

    return None



def build_isopy_obj(prob_values, formula_dict):

    prob_h_h = prob_values[0]
    prob_h_n = prob_values[1]
    prob_h_c = prob_values[2]
    prob_l_h = 1 - prob_h_h
    prob_l_n = 1 - prob_h_n
    prob_l_c = 1 - prob_h_c

    isopyobj = IsoSpecPy.IsoTotalProb(0.8, formula_dict['Light Molecular Formula'],
                                        atomCounts=(formula_dict['Heavy C'], formula_dict['Heavy H'], formula_dict['Heavy N']),
                                        isotopeProbabilities=([prob_l_c,prob_h_c],[prob_l_h,prob_h_h],[prob_l_n,prob_h_n]),
                                        isotopeMasses=(isotope_mass["C"],isotope_mass["H"],isotope_mass['N']))
    
    return isopyobj


def objective_function_wasserstein(initial_values, formula_dict, observed_masses, observed_intensities, mass_tolerance):
    # Generate theoretical isotope distribution
    isopyobj = build_isopy_obj(initial_values, formula_dict)
    #isopyobj = isopyobj.binned(0.05)
    isopyobj.normalize()
    theoretical_masses, theoretical_probs = isopyobj.np_masses(), isopyobj.np_probs()
    
    # Extract and normalize observed peaks
    extracted_masses, normalized_observed_intensities = extract_peptide_peaks(
        observed_masses, observed_intensities, theoretical_masses, mass_tolerance
    )
    
    # There was no peptide found in the spectra, return a high number
    if np.sum(normalized_observed_intensities) == 0:
        return 1000
    # Calculate the difference (e.g., sum of squared differences)
    distance = wasserstein_distance(theoretical_masses, extracted_masses, theoretical_probs, normalized_observed_intensities)

    return distance

def objective_function(initial_values, formula_dict, observed_masses, observed_intensities, mass_tolerance):
    # Generate theoretical isotope distribution
    isopyobj = build_isopy_obj(initial_values, formula_dict)
    #isopyobj = isopyobj.binned(0.05)
    isopyobj.normalize()
    theoretical_masses, theoretical_probs = isopyobj.np_masses(), isopyobj.np_probs()
    
    # Extract and normalize observed peaks
    extracted_masses, normalized_observed_intensities = extract_peptide_peaks(
        observed_masses, observed_intensities, theoretical_masses, mass_tolerance
    )
    
    theoretical_probs = theoretical_probs / np.max(theoretical_probs)
    
    # Calculate the difference (e.g., sum of squared differences)
    difference = np.sum((normalized_observed_intensities - theoretical_probs) ** 2)
    return difference

def plot_best_fit(best_values, formula_dict, observed_masses, observed_intensities, mass_tolerance):
    # Generate theoretical isotope distribution
    best_values = [best_values[i] for i in best_values]
    
    isopyobj = build_isopy_obj(best_values, formula_dict)
    #isopyobj = isopyobj.binned(0.05)
    isopyobj.normalize()
    
    theoretical_masses, theoretical_probs = isopyobj.np_masses(), isopyobj.np_probs()
    # Normalize to max intensity
    theoretical_probs = theoretical_probs / np.max(theoretical_probs)
    
    # Extract and normalize observed peaks
    extracted_masses, normalized_observed_intensities = extract_peptide_peaks(
        observed_masses, observed_intensities, theoretical_masses, mass_tolerance
    )
    compare_theo_obs_graphically(theoretical_masses, theoretical_probs, extracted_masses, normalized_observed_intensities, formula_dict['Peptide Sequence'])

def fit_isotope_abundances(formula_dict, observed_masses, observed_intensities, initial_guess):
    # Bounds for isotope abundances (between 0 and 1)
    bounds = [(0.9, 0.999), (0.97,0.999), (0.97,0.999)]
    # Convert initial_guess dict to a list
    initial_values = [initial_guess[element] for element in initial_guess]
    
    # Optimization
    result = minimize(
        objective_function_wasserstein,
        initial_values,
        args=(formula_dict, observed_masses, observed_intensities, 50),
        bounds=bounds,
        method='L-BFGS-B'
    )
    
    # Map the result back to isotope abundances
    optimized_abundances = dict(zip(initial_guess.keys(), result.x))
    return optimized_abundances

    
# Example sequence
sequence = "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"

# Perform the theoretical digestion
theoretical_df = theoretical_digest(sequence)
theoretical_df = theoretical_df.drop_duplicates(subset=['Peptide Sequence'])

# Find exchangeable protons
theoretical_df = determine_ex_H(theoretical_df)

# Decompose the formulas for isospecpy
theoretical_df = subtract_heavy_atoms(theoretical_df)

ebd = ElementDB()
hydrogen = ebd.getElement("H")
oxygen = ebd.getElement("O")
sulfur = ebd.getElement("S")
carbon = ebd.getElement("C")
nitrogen = ebd.getElement("N")

# Define natural abundances
isotope_dis = {
            "H":[i.getIntensity() for i in hydrogen.getIsotopeDistribution().getContainer()][0:2],
            "O":[i.getIntensity() for i in oxygen.getIsotopeDistribution().getContainer()],
            "S":[i.getIntensity() for i in sulfur.getIsotopeDistribution().getContainer()],
            "C":[i.getIntensity() for i in carbon.getIsotopeDistribution().getContainer()],
            "N":[i.getIntensity() for i in nitrogen.getIsotopeDistribution().getContainer()],
            }

# Define masses of isotopes
isotope_mass = {
            "H":[i.getMZ() for i in hydrogen.getIsotopeDistribution().getContainer()][0:2],
            "O":[i.getMZ() for i in oxygen.getIsotopeDistribution().getContainer()],
            "S":[i.getMZ() for i in sulfur.getIsotopeDistribution().getContainer()],
            "C":[i.getMZ() for i in carbon.getIsotopeDistribution().getContainer()],
            "N":[i.getMZ() for i in nitrogen.getIsotopeDistribution().getContainer()],
            }

exp_file = sys.argv[1]

initial_guess = {"H": 0.94, "N":0.98, "C":0.98}


exp_spectrum = GetPeakList(exp_file, 12)
obs_peak_values, obs_peak_intensities = exp_spectrum.get_peaks()

for row in theoretical_df.iterrows():
    peptide_dict = row[1].to_dict()
    optimized_abundances = fit_isotope_abundances(
        peptide_dict, obs_peak_values, obs_peak_intensities, initial_guess
    )

    print(f"Optimized Isotope Abundances for Peptide {peptide_dict['Peptide Sequence']}:", optimized_abundances)
    plot_best_fit(optimized_abundances, peptide_dict, obs_peak_values, obs_peak_intensities, 50)