""" Test for mdrun cli script

"""

import os
import subprocess

from aiida.orm.nodes.process.process import ProcessState

from aiida_gromacs.utils import searchprevious

from . import TEST_DIR


def test_launch_mdrun():
    """
    Run an instance of mdrun.
    """
    # get input file paths
    tpr_file = os.path.join(TEST_DIR, "input_files", "mdrun_1AKI_em.tpr")

    subprocess.check_output(
        [
            "gmx_mdrun",
            "-s",
            tpr_file,
            "-c",
            "1AKI_minimised.gro",
            "-e",
            "1AKI_minimised.edr",
            "-g",
            "1AKI_minimised.log",
            "-o",
            "1AKI_minimised.trr",
        ]
    )
    # append run process to qb
    # pylint: disable=unused-variable
    qb = searchprevious.build_query()
    # check the process has completed first
    # searchprevious.check_prev_process(qb)
    # pylint: disable=unsubscriptable-object
    prev_calc = qb.first()[0]

    if prev_calc.exit_status != 0:
        subprocess.check_output(
            [
                "gmx_mdrun",
                "-s",
                tpr_file,
                "-c",
                "1AKI_minimised.gro",
                "-e",
                "1AKI_minimised.edr",
                "-g",
                "1AKI_minimised.log",
                "-o",
                "1AKI_minimised.trr",
                "-ntomp",
                "5",
                "-ntmpi",
                "1",
            ]
        )

    # append run process to qb
    # pylint: disable=unused-variable
    qb = searchprevious.build_query()
    # check the process has completed first
    # searchprevious.check_prev_process(qb)
    # pylint: disable=unsubscriptable-object
    prev_calc = qb.first()[0]
    # check the process has finished and exited correctly
    assert prev_calc.process_state == ProcessState.FINISHED
    assert prev_calc.exit_status == 0
