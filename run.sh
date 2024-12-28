#!/bin/bash
set -xe

echo "Generating Topology"
curl -O https://files.rcsb.org/download/1AKI.pdb
grep -v HOH 1AKI.pdb > 1AKI_clean.pdb
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce <<EOF
15
EOF
echo "------------------------------------------------------------------------"

echo "Simulating a simple aqueous system"
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
echo "------------------------------------------------------------------------"

echo "Add ions..."
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral <<EOF
13
EOF
echo "------------------------------------------------------------------------"

echo "Energy Minimization..."
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
echo "------------------------------------------------------------------------"

echo "Equilibration - Phase 1..."
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
echo "------------------------------------------------------------------------"

echo "Equilibration - Phase 2..."
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt
echo "------------------------------------------------------------------------"

echo "Production MD"
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -v -deffnm md_0_1
echo "------------------------------------------------------------------------"


### Running GROMACS ON GPU
# gmx mdrun -deffnm md_0_1 -nb gp

## Analysis
# gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
# gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
# gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
#
