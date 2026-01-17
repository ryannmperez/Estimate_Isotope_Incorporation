# Isotope Distribution Estimation

Tools for calculating and visualizing fine isotope patterns in MALDI-TOF mass spectrometry data. Estimates heavy isotope incorporation fractions (13C, 15N, 2H) in tryptic peptides using Wasserstein distance optimization.

![Isotope Fitting Animation](output/EGVVHGVATVAEK_optimization.gif)

## Features

- **Spectrum Processing Pipeline**: Baseline subtraction, Savitzky-Golay smoothing, peak picking, and normalization
- **In-silico Digestion**: Supports Trypsin, GluC, and other proteases via pyopenms
- **Isotope Pattern Fitting**: Uses IsoSpecPy with threshold-based peak selection for accurate fitting across all peptide sizes
- **Optimization**: L-BFGS-B minimization of Wasserstein distance between theoretical and observed isotope distributions
- **Visualization**: Animated GIFs showing the optimization process with dark-themed styling
- **Batch Processing**: Parallel processing of multiple peptides

## Installation

```bash
pip install pyopenms IsoSpecPy numpy pandas matplotlib scipy biopython
```

## Usage

```bash
python calculate_fine_isotope_pattern.py <spectrum_file> --fasta <protein.fasta> [options]
```

### Basic Example

```bash
python calculate_fine_isotope_pattern.py spectrum.mzXML \
    --fasta protein.fasta \
    --enzyme Trypsin \
    --ppm 100 \
    --mode pos \
    --initial-h 94 \
    --initial-c 99 \
    --initial-n 99 \
    --animate \
    --save-stages \
    -o ./output
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `input_file` | Input spectrum file (.mzML or .mzXML) | Required |
| `--fasta, -f` | FASTA file containing protein sequence | Required |
| `--enzyme, -e` | Digestion enzyme (Trypsin, GluC, etc.) | Trypsin |
| `--ppm, -p` | PPM tolerance for peak matching | 50 |
| `--mode, -m` | Ionization mode (pos/neg) | pos |
| `--bin-size, -b` | Bin size for isotope distribution | 0.05 |
| `--initial-h` | Initial guess for 2H incorporation (%) | 94 |
| `--initial-c` | Initial guess for 13C incorporation (%) | 98.9 |
| `--initial-n` | Initial guess for 15N incorporation (%) | 98.9 |
| `--min-h` | Minimum bound for 2H optimization (%) | 90 |
| `--max-h` | Maximum bound for 2H optimization (%) | 99.9 |
| `--iso-threshold` | IsoSpecPy probability threshold | 0.001 |
| `--intensity-threshold` | Minimum relative peak intensity | 0.01 |
| `--mass-window` | Mass window around expected heavy mass (Da) | 15 |
| `--animate, -a` | Generate optimization animation GIFs | False |
| `--save-stages` | Save spectrum processing stage plots/mzML | False |
| `--output, -o` | Output directory | . |
| `--verbose, -v` | Verbose output with diagnostics | False |

## Output

- `isotope_incorporation_results.csv` - Results sorted by fit quality (Wasserstein distance)
- `*_optimization.gif` - Animation of optimization process (with `--animate`)
- `spectrum_*.png` - Processing stage plots (with `--save-stages`)
- `spectrum_*.mzML` - Processed spectra at each stage (with `--save-stages`)

## Methodology

### Spectrum Processing
1. **Baseline Subtraction**: Morphological filter to remove baseline
2. **Smoothing**: Savitzky-Golay filter for noise reduction
3. **Centroiding**: Iterative peak picking with S/N threshold
4. **Normalization**: Scale intensities to maximum of 1.0

### Isotope Pattern Calculation
Uses IsoSpecPy's `IsoThreshold` method to compute theoretical isotope distributions:
- Separates exchangeable H (natural abundance) from non-exchangeable H (can be deuterated)
- Models 13C, 15N, and 2H incorporation independently
- Threshold-based peak selection handles wide distributions from large peptides

### Optimization
- Minimizes Wasserstein (Earth Mover's) distance between theoretical and observed patterns
- L-BFGS-B bounded optimization
- Automatically filters observed peaks by intensity and mass window

## Example Results

```
TVEGAGSIAAATGFVK:
  Peaks in window: 34
  Wasserstein Distance: 0.2509
  Heavy H: 94.1%  Heavy N: 99.0%  Heavy C: 99.0%

EGVVHGVATVAEK:
  Peaks in window: 25
  Wasserstein Distance: 0.2635
  Heavy H: 95.1%  Heavy N: 99.2%  Heavy C: 99.2%
```

## License

MIT License
