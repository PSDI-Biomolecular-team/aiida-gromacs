#!/bin/bash

./launch.py --code gmx@localhost --command "pdb2gmx -i 1AKI_restraints.itp -o 1AKI_forcfield.gro -p 1AKI_topology.top -ff oplsaa -water spce -f 1AKI_clean.pdb" --inputs gromacs_files/1AKI_clean.pdb --outputs 1AKI_restraints.itp --outputs 1AKI_topology.top --outputs 1AKI_forcfield.gro


./launch.py --code gmx@localhost --command "editconf -f 1AKI_forcfield.gro -center 0 0 0 -d 1.0 -bt cubic -o 1AKI_newbox.gro" --inputs outputs/1AKI_forcfield.gro --outputs 1AKI_newbox.gro


./launch.py --code gmx@localhost --command "solvate -cp 1AKI_newbox.gro -p 1AKI_topology.top -cs spc216.gro -o 1AKI_solvated.gro" --inputs outputs/1AKI_newbox.gro --inputs outputs/1AKI_topology.top --outputs 1AKI_solvated.gro --outputs 1AKI_topology.top


./launch.py --code gmx@localhost --command "grompp -o 1AKI_ions.tpr -f ions.mdp -c 1AKI_solvated.gro -p 1AKI_topology.top" --inputs gromacs_files/ions.mdp --inputs outputs/1AKI_topology.top --inputs outputs/1AKI_solvated.gro --outputs 1AKI_ions.tpr
