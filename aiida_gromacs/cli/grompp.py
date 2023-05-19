#!/usr/bin/env python
"""CLI utility to run gmx grompp with AiiDA.

Usage: gmx_grompp --help
"""

import os

import click

from aiida import cmdline, engine
from aiida.plugins import CalculationFactory, DataFactory

from aiida_gromacs import helpers


def launch(params):
    """Run grompp.

    Uses helpers to add gromacs on localhost to AiiDA on the fly.
    """

    # If code is not initialised, then setup.
    gromacs_code = params.pop("code")
    if not gromacs_code:
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point="gromacs", computer=computer)

    # Prepare input parameters in AiiDA formats.
    SinglefileData = DataFactory("core.singlefile")
    mdpfile = SinglefileData(file=os.path.join(os.getcwd(), params.pop("f")))
    grofile = SinglefileData(file=os.path.join(os.getcwd(), params.pop("c")))
    topfile = SinglefileData(file=os.path.join(os.getcwd(), params.pop("p")))

    GromppParameters = DataFactory("gromacs.grompp")
    parameters = GromppParameters(params)

    # set up calculation
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "mdpfile": mdpfile,
        "grofile": grofile,
        "topfile": topfile,
        "metadata": {
            "description": "grompp job submission with the aiida_gromacs plugin",
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # pylint: disable=unused-variable
    future = engine.submit(CalculationFactory("gromacs.grompp"), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
@click.option("-f", default="grompp.mdp", type=str, help="Input parameter file")
@click.option("-c", required=True, type=str, help="Input structure file")
@click.option("-p", default="topol.top", type=str, help="Topology file")
@click.option("-o", default="conf.gro", type=str, help="Output structure file")
def cli(*args, **kwargs):
    # pylint: disable=unused-argument
    """Run example.

    Example usage:

    $ gmx_grompp --code gmx@localhost -f ions.mdp -c 1AKI_solvated.gro -p 1AKI_topology.top -o 1AKI_ions.tpr

    Alternative (automatically tried to create gmx@localhost code, but requires
    gromacs to be installed and available in your environment path):

    $ gmx_grompp -f ions.mdp -c 1AKI_solvated.gro -p 1AKI_topology.top -o 1AKI_ions.tpr

    Help: $ gmx_grompp --help
    """

    launch(kwargs)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
