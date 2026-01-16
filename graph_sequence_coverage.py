import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np
import sys

sequence = "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"

# Convert the sequence into a list of residues
sequence_residues = list(sequence)
sequence_length = len(sequence_residues)

# Get a list of CSV files from command-line arguments
file_names = sys.argv[1:]

if not file_names:
    print("Usage: python script_name.py <csv_file1> <csv_file2> ...")
    sys.exit(1)

# Initialize an array to hold residue coverage
coverage_array = np.zeros(sequence_length)

# Function to find all occurrences of a peptide in the full sequence
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
        # Look for the next occurrence
        index = full_sequence.find(peptide, index + 1)
    return positions

# Iterate over each CSV file provided
for file_name in file_names:
    # Read the CSV file
    df = pd.read_csv(file_name)
    # Extract unique peptides
    peptides = df['Peptide Sequence'].unique()
    # For each peptide, find positions and update coverage_array
    for peptide in peptides:
        positions = find_peptide_positions(peptide, sequence)
        if not positions:
            print(f"Warning: Peptide '{peptide}' not found in sequence.")
        for pos in positions:
            coverage_array[pos - 1] += 1  # Adjust for zero-based indexing

# Determine the maximum coverage value
max_coverage = int(np.max(coverage_array))

# Create a custom colormap
cmap_colors = ['grey', 'palegreen', 'forestgreen', 'darkgreen', 'navy', 'purple', 'red', 'yellow', 'orange', 'black']
# Select enough colors based on max_coverage
if max_coverage + 1 > len(cmap_colors):
    # Use a colormap with enough colors if max_coverage exceeds predefined colors
    cmap = plt.cm.get_cmap('viridis', max_coverage + 1)
else:
    cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', cmap_colors[:max_coverage + 1])

# Create a figure with reduced height
fig, ax = plt.subplots(figsize=(15, 2))

# Create a heatmap-like plot
data = coverage_array.reshape(1, -1)
im = ax.imshow(data, aspect='auto', cmap=cmap)

# Set major ticks at every 10 residues
major_xticks = np.arange(9, sequence_length, 10)
ax.set_xticks(major_xticks)
ax.set_xticklabels(np.arange(10, sequence_length + 1, 10), fontsize=10, fontweight='bold')

# Set minor ticks at every residue position
minor_xticks = np.arange(sequence_length)
ax.set_xticks(minor_xticks, minor=True)
ax.set_xticklabels(sequence_residues, fontsize=8, fontweight='bold', minor=True)

# Adjust tick parameters
ax.tick_params(axis='x', which='major', pad=15)  # Move numbers away from the axis
ax.tick_params(axis='x', which='minor', pad=5)   # Keep letters close to the axis

# Remove y-axis ticks
ax.set_yticks([])

# Adjust layout to make space for labels
plt.subplots_adjust(bottom=0.3)

# Add a colorbar with better formatting
cbar = fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.03, pad=0.04)
cbar.set_label('Residue Coverage', labelpad=20, fontsize=10, fontweight='bold')
# Set colorbar ticks and labels
cbar_ticks = np.arange(0, max_coverage + 1)
cbar.set_ticks(cbar_ticks)
cbar_labels = ['Not Covered'] + [f'Matches {i} Time' if i == 1 else f'Matches {i} Times' for i in cbar_ticks[1:]]
cbar.set_ticklabels(cbar_labels)
cbar.ax.tick_params(labelsize=9)

# Set title with enhanced style
ax.set_title('Peptide Matches Within ± 0.5% of NMR', fontsize=16, fontweight='bold')

# Adjust layout for better spacing
plt.tight_layout()

# Display the plot
plt.show()
