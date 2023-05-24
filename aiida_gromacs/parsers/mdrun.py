"""
Parsers provided by aiida_gromacs.

This calculation configures the ability to use the 'gmx mdrun' executable.
"""
from aiida.common import exceptions
from aiida.engine import ExitCode
from aiida.orm import SinglefileData
from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory

MdrunCalculation = CalculationFactory("gromacs.mdrun")


class MdrunParser(Parser):
    """
    Parser class for parsing output of calculation.
    """

    def __init__(self, node):
        """
        Initialize Parser instance

        Checks that the ProcessNode being passed was produced by a MdrunCalculation.

        :param node: ProcessNode of calculation
        :param type node: :class:`aiida.orm.nodes.process.process.ProcessNode`
        """
        super().__init__(node)
        if not issubclass(node.process_class, MdrunCalculation):
            raise exceptions.ParsingError("Can only parse MdrunCalculation")

    def parse(self, **kwargs):
        """
        Parse outputs, store results in database.

        :returns: an exit code, if parsing fails (or nothing if parsing succeeds)
        """
        outputs = ["stdout", "grofile", "enfile", "logfile", "trrfile"]

        # Check that folder content is as expected
        files_retrieved = self.retrieved.base.repository.list_object_names()
        files_expected = [
            self.node.get_option("output_filename"),
            self.node.inputs.parameters["c"],
            self.node.inputs.parameters["e"],
            self.node.inputs.parameters["g"],
            self.node.inputs.parameters["o"],
        ]

        if "cpo" in self.node.inputs.parameters.keys():
            outputs.append("cptfile")
            files_expected.append(self.node.inputs.parameters["cpo"])

        # Note: set(A) <= set(B) checks whether A is a subset of B
        if not set(files_expected) <= set(files_retrieved):
            self.logger.error(
                f"Found files '{files_retrieved}', expected to find '{files_expected}'"
            )
            return self.exit_codes.ERROR_MISSING_OUTPUT_FILES

        with self.retrieved.base.repository.open("_scheduler-stderr.txt", "rb") as handle:
             print(handle.read())

        # add outputs
        for index, thing in enumerate(files_expected):
            self.logger.info(f"Parsing '{thing}'")
            with self.retrieved.base.repository.open(thing, "rb") as handle:
                output_node = SinglefileData(filename=thing, file=handle)
            self.out(outputs[index], output_node)

        return ExitCode(0)
