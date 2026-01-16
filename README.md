# Isotope Distribution Analysis for Heavy-Labeled Peptides

Analysis tools for fitting isotope distributions of heavy-labeled (13C, 15N, 2H) peptides from mass spectrometry data, with comparison to NMR-derived deuteration fractions.

## Overview

This project provides tools to:
1. Theoretically digest protein sequences with different enzymes (Trypsin, GluC)
2. Calculate theoretical isotope patterns for heavy-labeled peptides
3. Fit observed mass spectra to determine isotope incorporation levels
4. Compare MS-derived deuteration with NMR measurements
5. Visualize sequence coverage and optimization results

## Requirements

```bash
pip install pyopenms IsoSpecPy numpy pandas matplotlib scipy
```

## Scripts

### `calculate_fine_isotope_pattern_simple.py`
Simple grid-search approach for fitting hydrogen deuteration from MS data. Compares results with NMR deuteration fractions.

```bash
python calculate_fine_isotope_pattern_simple.py <mzXML/mzML file> <enzyme> <ppm> <ion_mode> <bin_size>
```

**Arguments:**
- `mzXML/mzML file`: Input mass spectrometry file
- `enzyme`: `trypsin` or `gluc` (glutamyl endopeptidase)
- `ppm`: Mass tolerance in ppm (e.g., 30)
- `ion_mode`: `pos` or `neg`
- `bin_size`: Bin width for isotope pattern (e.g., 0.3)

**Example:**
```bash
python calculate_fine_isotope_pattern_simple.py sample.mzXML trypsin 30 pos 0.3
```

### `calculate_fine_isotope_pattern_with_animation.py`
Extended version that creates animated GIFs showing the optimization process for each peptide.

```bash
python calculate_fine_isotope_pattern_with_animation.py <mzXML/mzML file> <enzyme> <ppm> <ion_mode> <bin_size>
```

### `calculate_fine_isotope_pattern.py`
Original implementation using Wasserstein distance for optimization.

### `calculate_fine_isotope_pattern_with_sodium.py`
Handles both protonated [M+H]+ and sodiated [M+Na]+ adducts.

### `graph_sequence_coverage.py`
Visualizes peptide sequence coverage from processed CSV files.

```bash
python graph_sequence_coverage.py <processed_csv_1> <processed_csv_2> ...
```

## Input Data

### NMR Reference Data
- `Asyn_deuteration_fraction_CA.txt`: CA (alpha carbon) deuteration fractions from NMR
- `Asyn_deuteration_fraction_CB.txt`: CB (beta carbon) deuteration fractions from NMR

### Peptide Lists
- `Trypsin_Matching_Peptides.csv`: Expected peptides from trypsin digestion
- `GluC_Matching_Peptides.csv`: Expected peptides from GluC digestion

## Output

The scripts generate:
- `sorted_peptides.csv`: All peptides sorted by fit quality
- `<filename>_processed.csv`: Filtered results with MS vs NMR comparison
- `*_Theo_Exp.png`: Theoretical vs experimental isotope pattern plots
- `*_optimization_animation.gif`: Optimization process animations (with_animation version)

## Methodology

### Isotope Pattern Calculation
Uses IsoSpecPy to compute fine isotope patterns accounting for:
- 13C labeling (all carbons heavy-labeled)
- 15N labeling (all nitrogens heavy-labeled)
- 2H labeling (non-exchangeable hydrogens heavy-labeled)

Exchangeable hydrogens (backbone amides, side-chain OH/NH groups) are modeled separately.

### Optimization
Fits the heavy isotope abundance (primarily 2H) by minimizing the difference between theoretical and observed isotope patterns using:
- Grid search (simple version)
- L-BFGS-B optimization (animation version)
- Wasserstein distance metric

### NMR Comparison
Theoretical deuteration percentages are calculated from NMR CA/CB deuteration data, accounting for:
- Position-specific deuteration when available
- Amino acid-averaged deuteration for missing positions
- Side-chain deuteration estimates

## Target Protein

Default protein sequence is alpha-synuclein (140 residues):
```
MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA
```

Modify the `sequence` variable in the scripts to analyze different proteins.

## License

[Add your license here]

## Citation

[Add citation information if applicable]
