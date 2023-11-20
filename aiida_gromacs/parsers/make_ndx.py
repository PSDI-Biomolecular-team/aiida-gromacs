"""
Parsers provided by aiida_gromacs.

This calculation configures the ability to use the 'gmx make_ndx' executable.
"""
import os
from aiida.common import exceptions
from aiida.engine import ExitCode
from aiida.orm import SinglefileData
from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory

Make_ndxCalculation = CalculationFactory("gromacs.make_ndx")


class Make_ndxParser(Parser):
    """
    Parser class for parsing output of calculation.
    """

    def __init__(self, node):
        """
        Initialize Parser instance

        Checks that the ProcessNode being passed was produced by a Make_ndxCalculation.

        :param node: ProcessNode of calculation
        :param type node: :class:`aiida.orm.nodes.process.process.ProcessNode`
        """
        super().__init__(node)
        if not issubclass(node.process_class, Make_ndxCalculation):
            raise exceptions.ParsingError("Can only parse Make_ndxCalculation")

    def parse(self, **kwargs):
        """
        Parse outputs, store results in database.

        :returns: an exit code, if parsing fails (or nothing if parsing succeeds)
        """
        # Map output files to how they are named.
        outputs = ["stdout"]
        output_template = {
                "o": "n_file_out",
            }

        for item in output_template:
            if item in self.node.inputs.parameters.keys():
                outputs.append(output_template[item])

        # Grab list of retrieved files.
        files_retrieved = self.retrieved.base.repository.list_object_names()

        # Grab list of files expected and remove the scheduler stdout and stderr files.
        files_expected = [files for files in self.node.get_option("retrieve_list") if files not in ["_scheduler-stdout.txt", "_scheduler-stderr.txt"]]

        # Check if the expected files are a subset of retrieved.
        if not set(files_expected) <= set(files_retrieved):
            self.logger.error(
                f"Found files '{files_retrieved}', expected to find '{files_expected}'"
            )
            return self.exit_codes.ERROR_MISSING_OUTPUT_FILES

        # Map retrieved files to data nodes.
        for i, f in enumerate(files_expected):
            self.logger.info(f"Parsing '{f}'")
            with self.retrieved.base.repository.open(f, "rb") as handle:
                output_node = SinglefileData(filename=f, file=handle)
            self.out(outputs[i], output_node)

        # If not in testing mode, then copy back the files.
        if "PYTEST_CURRENT_TEST" not in os.environ:
            self.retrieved.copy_tree(os.getcwd())

        return ExitCode(0)