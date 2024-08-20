---

title: 'aiida-gromacs: A python plugin for dynamically tracking molecular dynamics simulations performed with GROMACS'
tags:
  - Python
  - data provenance
  - biomolecular simulation
  - reproducibility
  - FAIR data
  - computational biology
  - computation chemistry
authors:
  - name: Jas Kalayan
    orcid:  0000-0002-6833-1864
    equal-contrib: true
    affiliation: 1
  - name: James Gebbie-Rayet
    orcid:  0000-0001-8271-3431
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Science and Technology Facilities Council, Daresbury Laboratory, Sci-Tech Daresbury, UK
   index: 1
date: \today
bibliography: paper.bib

---

# Summary

[`aiida-gromacs`](https://aiida-gromacs.readthedocs.io/) is a python tool developed for capturing and sharing the full provenance for setup and performance of biomolecular simulations, tailored toward the popular MD engine [GROMACS](https://www.gromacs.org/) [@abraham_gromacs_2015]. It has been designed to mimic the way simulators currently work, allowing for rich metadata collection with minimal changes to how researchers do science.

# Statement of need

Popularity in molecular dynamics (MD) simulations has increased in research aimed at understanding biological functions and structural dynamics at atomistic scales. Improvements in forcefield potentials, faster performance from GPU processors and cheaper data storage enable the sampling of dynamics in larger and more complicated systems ($\sim 1 \times 10^5$ atoms producing $\sim$ GBs of data) at realistic compute timescales ($\sim$ hours).

Complementing computational research, experimental methods can now access increasing length scales from methods such as cryo- electron microscopy (EM) [@saibil_cryo-em_2022] and electron tomography (ET) [@neumuller_electron_2018] at near atomic resolutions. This has opened up research into using MD simulations to fit atomic structures into experimental density maps to assist in improving structure refinement tools [@joosten_roodmus_2024].

As more research communities require access to MD simulation data and how it has been curated, it is becoming increasingly important to appropriately share this data in a standardised format. Currently, there are limited tools available to track and collect workflows in molecular dynamics simulations of biomolecules. This is due to the complexity of workflows and the lack of standard practices to capture and disseminate them. The range of software used and the variety of simulation steps makes it difficult to keep track of biomolecular simulation protocols. Furthermore, simulation protocols are often difficult to reproduce from scientific publications, because typically the writing up of methods is carried out after the simulations have been performed. This means steps can be missing, mis-remembered or too convoluted to include in the resulting publication.

Of the limited tools available, there are several that aim to improve interoperability in simulation software using python scripts and libraries such as [BioBB](https://mmb.irbbarcelona.org/biobb/), which allows for simulations to be setup and run from external MD engines in a single python script [@andrio_bioexcel_2019], [MoSDeF](https://mosdef.org) is a collection of python libraries for building and parmeterising systems for simulations [@thompson_towards_2020] and [OpenMM](https://openmm.org) is used to setup and natively perform MD simulations within the python package [@eastman_openmm_2017].
Further to interoperability, [BioSimSpace](https://biosimspace.openbiosim.org/) allows for the creation of reusable components for setting up and performing simulations via [nodes](https://biosimspace.openbiosim.org/guides/nodes.html), which can be manually added for each step in a protocol to create a template workflow as a Jupyter notebook [@hedges_biosimspace_2019].

However, there are some limitations in current tools to produce truely [FAIR](https://www.go-fair.org/fair-principles/) [@wilkinson_fair_2016] biosimulation workflows, these include the assumed programming knowledge of researchers to create interoperable simulation workflows and the lack of standardisation of workflow formats for disseminating with the wider research community. We therefore use AiiDA, as it is an establised workflow tool that allows for the creation of plugins for any simulation software. Many of the current [registry](https://aiidateam.github.io/aiida-registry/) of AiiDA plugins are aimed at first-principles calculations that usually follow fixed workflows, yet, protocols in the computational biology field are tailored to the biological system of interest, with various software used to build and prepare molecular systems for simulations. Hence, `aiida-gromacs` has been specifically designed to flexibly capture a diverse range of simulation setup steps, allowing users to capture steps as they are performed, reducing the need to remember steps after the fact.

# Software overview

`aiida-gromacs` wraps native GROMACS commands to make it as easy as possible for the user to adapt their working, but are also allows for users to extend workflows to other programs that are run on the command line, meaning that no interoperability is required between software used. Full [documentation](https://aiida-gromacs.readthedocs.io),  along with several [tutorials](https://aiida-gromacs.readthedocs.io/en/latest/tutorials/) and Jupyter [notebooks](https://github.com/PSDI-UK/aiida-gromacs/tree/master/notebooks) for performing simulations via `aiida-gromacs` are available.

![A. Graph representation of provenance showing data nodes in green, commands submitted in red and the code used in blue. Connections between nodes are shown with arrows. B. Workflow for using aiida-gromacs and how it works once a command is submitted.\label{fig1:workflow}](figures/figure1.pdf)

## Command line interface

Using `aiida-gromacs` is as simple as submitting a command on your terminal, as all provenance is collected in the background via the AiiDA database. Commands native to GROMACS are wrapped into the `aiida-gromacs` plugin, the only difference in these commands is the added underscore between the `gmx` wrapper binary and the command, for example, `gmx mdrun` becomes `gmx_mdrun`. The same convention is used for other popular GROMACS commands run via the plugin. Additionally, we have included the `genericMD` tool to capture the provenance of commands used outside of GROMACS within a simulation workflow, allowing for simulators to produce unbroken provenance graphs of their complete simulation setup protocols.

Examples of how commands can be wrapped with `aiida-gromacs` are shown in the table below:

| Native command                   | `aiida-gromacs` wrapped command    |
| -------------------------------- | ---------------------------------- |
| `$ gmx mdrun -s mdrun.tpr -c mdrun.gro  -g mdrun.log -x mdrun.xtc` | `$ gmx_mdrun -s mdrun.tpr -c mdrun.gro  -g mdrun.log -x mdrun.xtc`|
| `$ martinize2 -f input.pdb -o out.top -x out.cg.pdb -ff martini3001 -nt -dssp mkdssp -elastic -p backbone`| `$ genericMD --code martinize2@localhost --command "-f input.pdb -o out.top -x out.cg.pdb -ff martini3001 -nt -dssp mkdssp -elastic -p backbone" --inputs input.pdb --outputs out.top --outputs out.cg.pdb --outputs molecule_0.itp`|

## Application programming interface

For simulators who prefer programming and want to include workflows in existing python scripts, provenance graphs can also be generated programmatically using the API.

### GROMACS specific commands
Import the relevent python libraries for using `aiida-core` and `aiida-gromacs`.
```python
from os import path
from aiida import engine, orm
from aiida.plugins import DataFactory, CalculationFactory
from aiida_gromacs import helpers
```

Store inputs as relevent datatypes such as the code, parameter and file datatypes. Each datatype is stored as a node in the provenance graph. `aiida-gromacs` handles the outputs automatically from included `gmx` commands and validates the correct outputs are produced.
```python
# 1. set up the GROMACS code datatype.
computer = helpers.get_computer()
gromacs_code = helpers.get_code(entry_point="gromacs",
                                computer=computer)

# 2. Prepare input parameters datatype
MdrunParameters = DataFactory("gromacs.mdrun")
parameters = MdrunParameters({"c": "mdrun.gro", "g": "mdrun.log",
                              "o": "mdrun.xtc", "v": "true"})

# 3. Define input files as AiiDA SinglefileData.
tprfile = orm.SinglefileData(file=path.join(os.getcwd(), "mdrun.tpr"))
```

Populate an input dictionary with all the data objects created above.
```python
# 4. Set up calculation dictionary
inputs = {
    "code": gromacs_code,
    "parameters": parameters,
    "tprfile": tprfile,
}
```

Run the job to run via `aiida-gromacs`, which functions as described in Figure \autoref{fig1:workflow}B.
```python
# 5. Run the calculation step in blocking mode.
result = engine.run(CalculationFactory("gromacs.mdrun"), **inputs)
```

### Generalise a CLI command with `genericMD`

As with CLI, the API methods for submitting jobs outside GROMACS are slightly modified with the `genericMD` tool.

Store inputs as relevent datatypes such as the code and file datatypes. Each datatype is stored as a node in the provenance graph.
```python
# 1. set up the GROMACS code datatype.
computer = helpers.get_computer()
gromacs_code = helpers.get_code(entry_point="martinize2",
                                computer=computer)


# 2. Define input files as AiiDA SinglefileData datatype.
input_files = {}
pdbfile = orm.SinglefileData(file=path.join(os.getcwd(), "input.pdb"))
input_files["pdbfile"] = pdbfile

# 3. Define output file names produced from the command,
# these must be known a priori
output_files = orm.List(["out.top", "out.cg.pdb", "molecule_0.itp"])

# 4. Save the command to run minus the code used to run it
command = orm.Str("-f input.pdb -o out.top -x out.cg.pdb "
            "-ff martini3001 -nt -dssp mkdssp -elastic -p backbone")
```

Populate an input dictionary with all the data objects created above.
```python
# 5. Set up calculation dictionary
inputs = {
    "code": gromacs_code,
    "command": command,
    "input_files": input_files,
    "output_files": output_files,
}
```

Run the job to run using the `genericMD` calculation factory in `aiida-gromacs`, which functions as described in Figure \autoref{fig1:workflow}B.
```python
# 6. Run the calculation step in blocking mode.
result = engine.run(CalculationFactory("genericMD"), **inputs)
```

## Outputs

The data provenance of simulation protocols are tracked and saved in a relational database. Inputs and outputs of commands can be represented as nodes in the graph and edges represent how data are connected (\autoref{fig1:workflow}).

Along with data provenance collection, the automation of simulation metadata collection is also included in `aiida-gromacs`. Output log files are automatically parsed and all relevant metadata associated with each simulation command is saved in json formatted dictionaries that can subsequently be used to automatically populate databases and vastly improve searchability of simulation data.


# Acknowledgements

PSDI, EPSRC, Archer2, HECBioSim, CPPBioSim

# References
