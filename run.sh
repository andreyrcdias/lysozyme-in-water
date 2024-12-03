#/bin/bash
set -x

## Generating Topology
# Step 1: Download the Protein Structure
echo "Downloading the Protein Structure 1AKI..."
curl -O https://files.rcsb.org/download/1AKI.pdb

# Step 2: Remove crystal clear water
echo "Removing Crystal Clear Water from the PDB file..."
grep -v HOH 1AKI.pdb > 1AKI_clean.pdb

# Step 3: Prepare the protein for simulation
echo "Preparing the protein for simulation..."
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce <<EOF
15
EOF
