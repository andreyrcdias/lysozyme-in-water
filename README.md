# lysozyme-in-water
[GROMACS Tutorial - Lysozyme in Water](http://www.mdtutorials.com/gmx/lysozyme/index.html)


## Prerequisites
- [Gromacs](https://www.gromacs.org/)
- [curl](https://curl.se/)


## Quick Start
> [!TIP]
> You can run the commands sequentially by your own or run the shell script:
```bash
./run.sh
```
> [!TIP]
> Also, if you want to clean the generated files, run:
```bash
./clean.sh
```
___


### Generate Topology
```bash
curl -O https://files.rcsb.org/download/1AKI.pdb
```
```bash
grep -v HOH 1AKI.pdb > 1AKI_clean.pdb
```
```bash
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce
```
> We will use the all-atom **`OPLS`** force field, so type **`15`** at the command prompt, followed by `Enter`


### Simulating a simple aqueous system
```bash
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
```
```bash
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
```


### Add Ions
```bash
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
```
```bash
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
```
```bash
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```
> When prompted, choose group **`13 "SOL"`** for embedding ions. You do not want to replace parts of your protein with ions.


### Energy Minimization
```bash
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
```
```bash
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
```
```bash
gmx mdrun -v -deffnm em
```


### Equilibration
#### Phase 1
```bash
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
```
```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```
```bash
gmx mdrun -v -deffnm nvt
```
#### Phase 2
```bash
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
```
```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
```
```bash
gmx mdrun -v -deffnm npt
```


### Production MD
```bash
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp
```
```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
```
```bash
gmx mdrun -v -deffnm md_0_1
```


### Analysis
```bash
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
```
```bash
gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
```
```bash
gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
```
