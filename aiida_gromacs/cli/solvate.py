#!/usr/bin/env python
"""CLI utility to run gmx solvate with AiiDA.

Usage: gmx_solvate --help
"""

import os

import click

from aiida import cmdline, engine
from aiida.plugins import CalculationFactory, DataFactory

from aiida_gromacs import helpers


def launch(params):
    """Run solvate.

    Uses helpers to add gromacs on localhost to AiiDA on the fly.
    """

    # If code is not initialised, then setup.
    gromacs_code = params.pop("code")
    if not gromacs_code:
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point="gromacs", computer=computer)

    # Prepare input parameters in AiiDA formats.
    SinglefileData = DataFactory("core.singlefile")
    grofile = SinglefileData(file=os.path.join(os.getcwd(), params.pop("cp")))
    topfile = SinglefileData(file=os.path.join(os.getcwd(), params.pop("p")))

    SolvateParameters = DataFactory("gromacs.solvate")
    parameters = SolvateParameters(params)

    # set up calculation
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "grofile": grofile,
        "topfile": topfile,
        "metadata": {
            "description": "solvate job submission with the aiida_gromacs plugin",
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # pylint: disable=unused-variable
    future = engine.submit(CalculationFactory("gromacs.solvate"), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
@click.option("-cp", default="protein.gro", type=str, help="Input structure file")
@click.option("-cs", default="spc216.gro", type=str, help="Library structure file")
@click.option("-p", default="topol.top", type=str, help="Topology file")
@click.option("-o", default="out.gro", type=str, help="Output structure file")
def cli(*args, **kwargs):
    # pylint: disable=unused-argument
    """Run example.

    Example usage:

    $ gmx_solvate --code gmx@localhost -cp 1AKI_newbox.gro -cs spc216.gro -p 1AKI_topology.top -o 1AKI_solvated.gro

    Alternative (automatically tried to create gmx@localhost code, but requires
    gromacs to be installed and available in your environment path):

    $ gmx_solvate -cp 1AKI_newbox.gro -cs spc216.gro -p 1AKI_topology.top -o 1AKI_solvated.gro

    Help: $ gmx_solvate --help
    """

    launch(kwargs)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
