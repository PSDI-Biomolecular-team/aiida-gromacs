#!/usr/bin/env python
"""CLI utility to run gmx mdrun with AiiDA.

Usage: gmx_mdrun --help
"""

import os

import click

from aiida import cmdline, engine
from aiida.plugins import CalculationFactory, DataFactory

from aiida_gromacs import helpers


def launch(params):
    """Run mdrun.

    Uses helpers to add gromacs on localhost to AiiDA on the fly.
    """

    # If code is not initialised, then setup.
    gromacs_code = params.pop("code")
    if not gromacs_code:
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point="gromacs", computer=computer)

    # Prepare input parameters in AiiDA formats.
    SinglefileData = DataFactory("core.singlefile")
    tprfile = SinglefileData(file=os.path.join(os.getcwd(), params.pop("s")))

    MdrunParameters = DataFactory("gromacs.mdrun")
    parameters = MdrunParameters(params)

    # set up calculation.
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "tprfile": tprfile,
        "metadata": {
            "description": "mdrun minimisation job submission with the aiida_gromacs plugin",
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # pylint: disable=unused-variable
    future = engine.submit(CalculationFactory("gromacs.mdrun"), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
@click.option("-s", default="topol.tpr", type=str, help="Portable xdr run input file")
@click.option("-c", default="confout.gro", type=str, help="Structure file")
@click.option("-e", default="ener.edr", type=str, help="Energy file")
@click.option("-g", default="md.log", type=str, help="MD log file")
@click.option("-o", default="topol.gro", type=str, help="Trajectory output file")
@click.option("-v", default="false", type=str, help="verbose")
def cli(*args, **kwargs):
    # pylint: disable=unused-argument
    # pylint: disable=line-too-long
    """Run example.

    Example usage:

    $ gmx_mdrun --code gmx@localhost -s 1AKI_em.tpr -c 1AKI_minimised.gro -e 1AKI_minimised.edr -g 1AKI_minimised.log -o 1AKI_minimised.trr

    Alternative (automatically tried to create gmx@localhost code, but requires
    gromacs to be installed and available in your environment path):

    $ gmx_mdrun -s 1AKI_em.tpr -c 1AKI_minimised.gro -e 1AKI_minimised.edr -g 1AKI_minimised.log -o 1AKI_minimised.trr

    Help: $ gmx_mdrun --help
    """

    launch(kwargs)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
