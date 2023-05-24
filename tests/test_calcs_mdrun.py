""" Tests for calculations

"""
import os

from aiida.engine import run
from aiida.plugins import CalculationFactory, DataFactory

from . import TEST_DIR


def test_process(gromacs_code):
    """Test running a mdrun calculation.
    Note: this does not test that the expected outputs are created of output parsing"""

    # Prepare input parameters
    MdrunParameters = DataFactory("gromacs.mdrun")
    parameters = MdrunParameters(
        {
            "c": "mdrun_1AKI_minimised.gro",
            "e": "mdrun_1AKI_minimised.edr",
            "g": "mdrun_1AKI_minimised.log",
            "o": "mdrun_1AKI_minimised.trr",
            "v": "true",
            "ntomp": "5",
            "ntmpi": "1",
        }
    )

    SinglefileData = DataFactory("core.singlefile")
    tprfile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "mdrun_1AKI_em.tpr")
    )

    # set up calculation
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "tprfile": tprfile,
        "metadata": {
            "description": "mdrun test",
        },
    }

    result = run(CalculationFactory("gromacs.mdrun"), **inputs)

    assert "stdout" in result
    assert "trrfile" in result
    assert "grofile" in result
    assert "logfile" in result
    assert "enfile" in result
