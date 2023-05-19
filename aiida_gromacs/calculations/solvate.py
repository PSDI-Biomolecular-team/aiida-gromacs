"""
Calculations provided by aiida_gromacs.

This calculation configures the ability to use the 'gmx solvate' executable.
"""
from aiida.common import CalcInfo, CodeInfo
from aiida.engine import CalcJob
from aiida.orm import SinglefileData
from aiida.plugins import DataFactory

SolvateParameters = DataFactory("gromacs.solvate")


class SolvateCalculation(CalcJob):
    """
    AiiDA calculation plugin wrapping the 'gmx solvate' executable.

    AiiDA plugin wrapper for solvating a molecular system.
    """

    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        # yapf: disable
        super().define(spec)

        # set default values for AiiDA options
        spec.inputs['metadata']['options']['resources'].default = {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 1,
        }
        spec.inputs['metadata']['options']['parser_name'].default = 'gromacs.solvate'
        spec.input('metadata.options.output_filename', valid_type=str, default='solvate.out')
        spec.input('grofile', valid_type=SinglefileData, help='Input structure')
        spec.input('topfile', valid_type=SinglefileData, help='Input topology')
        spec.input('parameters', valid_type=SolvateParameters, help='Command line parameters for gmx solvate.')

        spec.output('stdout', valid_type=SinglefileData, help='stdout')
        spec.output('grofile', valid_type=SinglefileData, help='Output solvated gro file.')
        spec.output('topfile', valid_type=SinglefileData, help='Output topology file.')

        spec.exit_code(300, 'ERROR_MISSING_OUTPUT_FILES', message='Calculation did not produce all expected output files.')

    def prepare_for_submission(self, folder):
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin should temporarily place all files
            needed by the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        codeinfo = CodeInfo()
        codeinfo.cmdline_params = self.inputs.parameters.cmdline_params(
            grofile=self.inputs.grofile.filename, topfile=self.inputs.topfile.filename
        )
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.withmpi = self.inputs.metadata.options.withmpi

        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = [
            (
                self.inputs.grofile.uuid,
                self.inputs.grofile.filename,
                self.inputs.grofile.filename,
            ),
            (
                self.inputs.topfile.uuid,
                self.inputs.topfile.filename,
                self.inputs.topfile.filename,
            ),
        ]
        calcinfo.retrieve_list = [
            self.metadata.options.output_filename,
            self.inputs.parameters["o"],
            self.inputs.topfile.filename,
        ]

        return calcinfo
