#!/bin/bash

launch --code gmx@localhost --command "pdb2gmx -i 1AKI_restraints.itp -o 1AKI_forcfield.gro -p 1AKI_topology.top -ff oplsaa -water spce -f 1AKI_clean.pdb" --inputs gromacs_files/1AKI_clean.pdb --outputs 1AKI_restraints.itp --outputs 1AKI_topology.top --outputs 1AKI_forcfield.gro


launch --code gmx@localhost --command "editconf -f 1AKI_forcfield.gro -center 0 0 0 -d 1.0 -bt cubic -o 1AKI_newbox.gro" --inputs outputs/1AKI_forcfield.gro --outputs 1AKI_newbox.gro


launch --code gmx@localhost --command "solvate -cp 1AKI_newbox.gro -p 1AKI_topology.top -cs spc216.gro -o 1AKI_solvated.gro" --inputs outputs/1AKI_newbox.gro --inputs outputs/1AKI_topology.top --outputs 1AKI_solvated.gro --outputs 1AKI_topology.top


launch --code gmx@localhost --command "grompp -o 1AKI_ions.tpr -f ions.mdp -c 1AKI_solvated.gro -p 1AKI_topology.top" --inputs gromacs_files/ions.mdp --inputs outputs/1AKI_topology.top --inputs outputs/1AKI_solvated.gro --outputs 1AKI_ions.tpr

launch --code bash@localhost --command "echo SOL | gmx genion -s 1AKI_ions.tpr  -p 1AKI_topology.top -o 1AKI_solvated_ions.gro -pname NA -nname CL -neutral true" --inputs outputs/1AKI_ions.tpr --inputs outputs/1AKI_topology.top --outputs 1AKI_solvated_ions.gro  --outputs 1AKI_topology.top

launch --code gmx@localhost --command "grompp -f min.mdp -c 1AKI_solvated_ions.gro -p 1AKI_topology.top -o 1AKI_min.tpr" --inputs gromacs_files/min.mdp --inputs outputs/1AKI_solvated_ions.gro --inputs outputs/1AKI_topology.top --outputs 1AKI_min.tpr

launch --code gmx@localhost --command "mdrun -s 1AKI_min.tpr -c 1AKI_minimised.gro -e 1AKI_minimised.edr -g 1AKI_minimised.log -o 1AKI_minimised.trr -v true " --inputs outputs/1AKI_min.tpr --outputs 1AKI_minimised.gro --outputs 1AKI_minimised.edr --outputs 1AKI_minimised.log --outputs 1AKI_minimised.trr
