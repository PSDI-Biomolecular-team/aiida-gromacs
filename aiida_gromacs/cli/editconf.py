#!/usr/bin/env python
"""CLI utility to run gmx editconf with AiiDA.

Usage: gmx_editconf --help
"""

import os

import click

from aiida import cmdline, engine
from aiida.plugins import CalculationFactory, DataFactory

from aiida_gromacs import helpers


def launch(params):
    """Run editconf.

    Uses helpers to add gromacs on localhost to AiiDA on the fly.
    """

    # If code is not initialised, then setup.
    gromacs_code = params.pop("code")
    if not gromacs_code:
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point="gromacs", computer=computer)

    # Prepare input parameters in AiiDA formats.
    SinglefileData = DataFactory("core.singlefile")
    grofile = SinglefileData(file=os.path.join(os.getcwd(), params.pop("f")))

    EditconfParameters = DataFactory("gromacs.editconf")
    parameters = EditconfParameters(params)

    # set up calculation.
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "grofile": grofile,
        "metadata": {
            "description": "record editconf data provenance via the aiida_gromacs plugin",
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # pylint: disable=unused-variable
    future = engine.submit(CalculationFactory("gromacs.editconf"), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
@click.option("-f", default="conf.gro", type=str, help="Input structure file")
@click.option("-center", default="0 0 0", type=str, help="Shift geometrical center")
@click.option("-d", default="0", type=str, help="Distance between box and solute")
@click.option("-bt", default="triclinic", type=str, help="Box type")
@click.option("-o", default="out.gro", type=str, help="Output structure file")
def cli(*args, **kwargs):
    # pylint: disable=unused-argument
    """Run example.

    Example usage:

    $ gmx_editconf --code gmx@localhost -f 1AKI_forcefield.gro -center 0 -d 1.0 -bt cubic -o 1AKI_newbox.gro

    Alternative (automatically tried to create gmx@localhost code, but requires
    gromacs to be installed and available in your environment path):

    $ gmx_editconf -f 1AKI_forcefield.gro -center 0 -d 1.0 -bt cubic -o 1AKI_newbox.gro

    Help: $ gmx_editconf --help
    """

    launch(kwargs)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
