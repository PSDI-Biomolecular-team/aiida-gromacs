""" Tests for calculations

"""
import os

from aiida.engine import run
from aiida.plugins import CalculationFactory, DataFactory

from . import TEST_DIR


def test_process(gromacs_code):
    """Test running a pdb2gmx calculation.
    Note: this does not test that the expected outputs are created of output parsing"""

    # Prepare input parameters
    Pdb2gmxParameters = DataFactory("gromacs.pdb2gmx")
    parameters = Pdb2gmxParameters(
        {
            "ff": "oplsaa",
            "water": "spce",
            "o": "pdb2gmx_1AKI_forcefield.gro",
            "p": "pdb2gmx_1AKI_topology.top",
            "i": "pdb2gmx_1AKI_restraints.itp",
        }
    )

    SinglefileData = DataFactory("core.singlefile")
    pdbfile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "pdb2gmx_1AKI_clean.pdb")
    )

    # set up calculation
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "pdbfile": pdbfile,
        "metadata": {
            "description": "pdb2gmx test",
        },
    }

    result = run(CalculationFactory("gromacs.pdb2gmx"), **inputs)

    import time
    time.delay(60)

    assert "stdout" in result
    assert "grofile" in result
    assert "itpfile" in result
    assert "topfile" in result
