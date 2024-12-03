# lysozyme-in-water

## Generate Topology
1. Download the Protein Structure [1AKI](https://www.rcsb.org/structure/1aki)
```bash
curl -O https://files.rcsb.org/download/1AKI.pdb
```
2. Removing crystal clear water with `grep`
```bash
grep -v HOH 1AKI.pdb > 1AKI_clean.pdb
```
3. Prepares the protein for simulation by creating essential files that define its structure and properties.
```bash
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce
```
> We will use the all-atom **`OPLS`** force field, so type **`15`** at the command prompt, followed by `Enter`


## Simulating a simple aqueous system

1. Defining box
```bash
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
```

2. Filling with solvent (water)
```bash
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
```


## Add Ions
1. Download the [IONS.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp)
```bash
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
```

2. Assemble your Ions `.tpr` file
```bash
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
```

3. Now we have an atomic-level description of our system in the binary file ions.tpr. We will pass this file to genion:
```bash
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```
> When prompted, choose group **`13 "SOL"`** for embedding ions. You do not want to replace parts of your protein with ions.


## Energy Minimization
1. Download the Energy Minimization (**EM**) [minim.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp)
```bash
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
```

2. Assemble the binary input using grompp using this input parameter file:
```bash
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
```

3. Invoke mdrun to carry out the EM:
```bash
gmx mdrun -v -deffnm em
```

<!-- Review if its needed -->
<!-- 4. Let's do a bit of analysis. The em.edr file contains all of the energy terms that GROMACS collects during EM. You can analyze any .edr file using the GROMACS energy module:

```bash
gmx energy -f em.edr -o potential.xvg
``` -->


## Equilibration

1. Download OPLS [nvt.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp)
```bash
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
```

2. We will call grompp and mdrun just as we did at the EM step:
```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```

<!-- Review if its needed -->
<!--
```bash
gmx mdrun -deffnm nvt
```
-->
<!-- 4. Analyze the temperature progression, again using energy: -->
<!--
```bash
gmx energy -f nvt.edr -o temperature.xvg
```
-->

3. Download OPLS [npt.mpd](http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp)
```bash
curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
```

___
[GROMACS Tutorial - Lysozyme in Water](http://www.mdtutorials.com/gmx/lysozyme/index.html)
