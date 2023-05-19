"""
Calculations provided by aiida_gromacs.

This calculation configures the ability to use the 'gmx mdrun' executable.
"""
from aiida.common import CalcInfo, CodeInfo
from aiida.engine import CalcJob
from aiida.orm import SinglefileData
from aiida.plugins import DataFactory

MdrunParameters = DataFactory("gromacs.mdrun")


class MdrunCalculation(CalcJob):
    """
    AiiDA calculation plugin wrapping the 'gmx mdrun' executable.

    AiiDA plugin wrapper for converting PDB files to GRO files.
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
            'num_cores_per_mpiproc': 5,
        }

        spec.inputs['metadata']['options']['max_wallclock_seconds'].default = 86400

        spec.inputs['metadata']['options']['parser_name'].default = 'gromacs.mdrun'
        spec.input('metadata.options.output_filename', valid_type=str, default='mdrun.out')
        spec.input('tprfile', valid_type=SinglefileData, help='Input structure.')
        spec.input('parameters', valid_type=MdrunParameters, help='Command line parameters for gmx mdrun')

        spec.output('stdout', valid_type=SinglefileData, help='stdout')
        spec.output('trrfile', valid_type=SinglefileData, help='Output trajectory.')
        spec.output('grofile', valid_type=SinglefileData, help='Output structure file.')
        spec.output('logfile', valid_type=SinglefileData, help='Output log file.')
        spec.output('enfile', valid_type=SinglefileData, help='Output energy file.')

        spec.output('cptfile', valid_type=SinglefileData, required=False, help='Checkpoint file.')

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
            tprfile=self.inputs.tprfile.filename
        )
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.withmpi = self.inputs.metadata.options.withmpi

        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = [
            (
                self.inputs.tprfile.uuid,
                self.inputs.tprfile.filename,
                self.inputs.tprfile.filename,
            ),
        ]

        calcinfo.retrieve_list = [
            self.metadata.options.output_filename,
            self.inputs.parameters["c"],
            self.inputs.parameters["e"],
            self.inputs.parameters["g"],
            self.inputs.parameters["o"],
        ]

        if "cpo" in self.inputs.parameters.keys():
            calcinfo.retrieve_list.append(self.inputs.parameters["cpo"])

        return calcinfo
