import pandas as pd
from pyopenms import *
import sys
import IsoSpecPy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Define the function for theoretical digestion
def theoretical_digest(sequence: str, enzyme: str = "Trypsin", missed_cleavages: int = 0):
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
    exp = MSExperiment()
    MzMLFile().load(file, exp)
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
    from pyopenms import EmpiricalFormula

    light_formulas = []
    heavy_C_list = []
    heavy_N_list = []
    heavy_H_list = []
    sodiated_formulas = []
    sodiated_heavy_H_list = []

    for idx, row in theoretical_df.iterrows():
        formula_str = row['Molecular Formula']
        ex_H = row['Ex. Hydrogens']

        formula = EmpiricalFormula(formula_str)
        element_counts = {}
        for elem in formula.getElementalComposition():
            symbol = elem.decode()
            count = formula.getElementalComposition()[elem]
            element_counts[symbol] = count

        H_total = element_counts.get('H', 0)
        H_heavy = H_total - ex_H

        C_heavy = element_counts.get('C', 0)
        N_heavy = element_counts.get('N', 0)

        heavy_C_list.append(C_heavy)
        heavy_N_list.append(N_heavy)
        heavy_H_list.append(H_heavy)

        H_light = ex_H
        new_element_counts = {}
        for elem, count in element_counts.items():
            if elem == 'C' or elem == 'N':
                continue
            elif elem == 'H':
                if H_light > 0:
                    new_element_counts['H'] = H_light
            else:
                new_element_counts[elem] = count

        new_formula = EmpiricalFormula()
        for elem, count in new_element_counts.items():
            new_formula += EmpiricalFormula(f"{elem}{int(count)}")

        light_formulas.append(new_formula.toString())

        # Generate sodiated formula
        sodiated_formula = EmpiricalFormula(new_formula.toString())
        sodiated_formula -= EmpiricalFormula("H1")
        sodiated_formula += EmpiricalFormula("Na1")
        sodiated_formulas.append(sodiated_formula.toString())

        # Adjust heavy hydrogens for sodiated form
        sodiated_heavy_H = H_heavy - 1
        sodiated_heavy_H_list.append(sodiated_heavy_H)

    theoretical_df['Light Molecular Formula'] = light_formulas
    theoretical_df['Heavy C'] = heavy_C_list
    theoretical_df['Heavy H'] = heavy_H_list
    theoretical_df['Heavy N'] = heavy_N_list
    theoretical_df['Sodiated Molecular Formula'] = sodiated_formulas
    theoretical_df['Sodiated Heavy H'] = sodiated_heavy_H_list
    theoretical_df['Sodiated Heavy C'] = heavy_C_list
    theoretical_df['Sodiated Heavy N'] = heavy_N_list

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
    else:
        # If total intensity is zero, return the original zero intensities
        normalized_intensities = extracted_intensities
    
    return extracted_masses, normalized_intensities

def compare_theo_obs_graphically(theo_masses, theo_probs, obs_masses, obs_probs, peptide):

    plt.bar(theo_masses, theo_probs, width=0.05, label='Theoretical Masses')
    plt.bar(obs_masses, obs_probs, width=0.05, label='Observed Masses')
    plt.legend()
    plt.title(f"{peptide} Theo. vs. Exp.")
    plt.xlabel('m/z')
    plt.ylabel("Normalized Intensity")
    plt.savefig(f"{peptide}_Theo_Exp.png")
    plt.close()
    
    return None


def build_isopy_obj(prob_values, formula_string, heavy_counts):
    formula = EmpiricalFormula(formula_string)
    element_counts = formula.getElementalComposition()

    elements = []
    atom_counts = []
    heavy_atom_counts = []
    isotope_probabilities = []
    isotope_masses = []

    for elem in sorted(element_counts.keys()):
        symbol = elem.decode()
        total_count = element_counts[elem]
        heavy_count = heavy_counts.get(symbol, 0)

        elements.append(symbol)
        atom_counts.append(total_count)
        heavy_atom_counts.append(heavy_count)

        masses = isotope_mass[symbol]
        probs = isotope_dis[symbol]

        if symbol in prob_values:
            prob_heavy = prob_values[symbol]
            prob_light = 1 - prob_heavy
            adjusted_probs = [prob_light, prob_heavy]
            adjusted_masses = masses[:2]
        else:
            adjusted_probs = probs
            adjusted_masses = masses

        isotope_probabilities.append(adjusted_probs)
        isotope_masses.append(adjusted_masses)

    isopyobj = IsoSpecPy.IsoTotalProb(
        0.9999,
        formula_string,
        atomCounts=heavy_atom_counts,
        isotopeProbabilities=isotope_probabilities,
        isotopeMasses=isotope_masses
    )
    return isopyobj




def objective_function(initial_values, elements, formula_string, heavy_counts, observed_masses, observed_intensities, mass_tolerance):
    # Construct prob_values dict from elements and initial_values
    prob_values = dict(zip(elements, initial_values))

    # Generate theoretical isotope distribution
    isopyobj = build_isopy_obj(prob_values, formula_string, heavy_counts)
    isopyobj.normalize()
    theoretical_masses, theoretical_probs = isopyobj.np_masses(), isopyobj.np_probs()

    # Extract and normalize observed peaks
    extracted_masses, normalized_observed_intensities = extract_peptide_peaks(
        observed_masses, observed_intensities, theoretical_masses, mass_tolerance
    )

    # Calculate the difference (e.g., sum of squared differences)
    difference = np.sum((normalized_observed_intensities - theoretical_probs) ** 2)
    return difference


def plot_best_fit(best_values, formula_dict, observed_masses, observed_intensities, mass_tolerance):
    # Generate theoretical isotope distribution
    best_values = [best_values[i] for i in best_values]
    isopyobj = build_isopy_obj(best_values, formula_dict)
    isopyobj.normalize()
    theoretical_masses, theoretical_probs = isopyobj.np_masses(), isopyobj.np_probs()
    
    # Extract and normalize observed peaks
    extracted_masses, normalized_observed_intensities = extract_peptide_peaks(
        observed_masses, observed_intensities, theoretical_masses, mass_tolerance
    )
    compare_theo_obs_graphically(theoretical_masses, theoretical_probs, extracted_masses, normalized_observed_intensities, formula_dict['Peptide Sequence'])

def fit_isotope_abundances(formula_string, heavy_counts, observed_masses, observed_intensities, initial_guess):
    initial_values_list = []
    bounds = []
    elements = []
    for elem, initial_value in initial_guess.items():
        elements.append(elem)
        initial_values_list.append(initial_value)
        bounds.append((0.9, 0.999))

    # Optimization
    result = minimize(
        objective_function,
        initial_values_list,
        args=(elements, formula_string, heavy_counts, observed_masses, observed_intensities, 50),
        bounds=bounds,
        method='L-BFGS-B'
    )

    # Map the result back to isotope abundances
    optimized_abundances = dict(zip(elements, result.x))
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
sodium = ebd.getElement("Na")



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

isotope_dis["Na"] = [i.getIntensity() for i in sodium.getIsotopeDistribution().getContainer()]
isotope_mass["Na"] = [i.getMZ() for i in sodium.getIsotopeDistribution().getContainer()]
exp_file = sys.argv[1]

initial_guess = {"H": 0.93, "N":0.99, "C":0.99}

exp_spectrum = GetPeakList(exp_file, 3)
obs_peak_values, obs_peak_intensities = exp_spectrum.get_peaks()

for idx, row in theoretical_df.iterrows():
    peptide_dict = row.to_dict()

    # Protonated form
    formula_string = peptide_dict['Light Molecular Formula']
    heavy_counts = {'C': peptide_dict['Heavy C'], 'H': peptide_dict['Heavy H'], 'N': peptide_dict['Heavy N']}
    optimized_abundances = fit_isotope_abundances(
        formula_string, heavy_counts, obs_peak_values, obs_peak_intensities, initial_guess
    )
    print(f"Optimized Isotope Abundances for Protonated Peptide {peptide_dict['Peptide Sequence']}:", optimized_abundances)
    #plot_best_fit(optimized_abundances, formula_string, heavy_counts, obs_peak_values, obs_peak_intensities, 50, peptide_dict['Peptide Sequence'] + "_protonated")

    # Sodiated form
    formula_string = peptide_dict['Sodiated Molecular Formula']
    heavy_counts = {'C': peptide_dict['Sodiated Heavy C'], 'H': peptide_dict['Sodiated Heavy H'], 'N': peptide_dict['Sodiated Heavy N']}
    optimized_abundances = fit_isotope_abundances(
        formula_string, heavy_counts, obs_peak_values, obs_peak_intensities, initial_guess
    )
    print(f"Optimized Isotope Abundances for Sodiated Peptide {peptide_dict['Peptide Sequence']}:", optimized_abundances)
    #plot_best_fit(optimized_abundances, formula_string, heavy_counts, obs_peak_values, obs_peak_intensities, 50, peptide_dict['Peptide Sequence'] + "_sodiated")
