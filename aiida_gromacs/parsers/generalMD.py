from aiida.engine import ExitCode
from aiida.orm import SinglefileData
from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory

import re
import os

# entry point string under which the parser class is registered:
GeneralCalculation = CalculationFactory('general-MD') 


class GeneralParser(Parser):
    """
    Parsing the output files produced by a code into AiiDA nodes 
        is optional, but it can make your data queryable and therefore easier 
        to access and analyze.
    Before the parse() method is called, two important attributes are set on 
        the Parser instance:
        1. self.retrieved: An instance of FolderData, which points to the 
            folder containing all output files that the CalcJob instructed 
            to retrieve, and provides the means to open() any file it contains.
        2. self.node: The CalcJobNode representing the finished calculation, 
            which, among other things, provides access to all of its inputs 
            (self.node.inputs).

    """
    def parse(self, **kwargs):
        """
        Parse outputs, store results in database.
        """

        # dirpath = pathlib.Path(kwargs['retrieved_temporary_folder'])
        # remote_folder = self.node.outputs["remote_folder"] #.listdir(path)
        
        # get_option() convenience method is used to get the filename of 
        # the output file.
        output_filename = self.node.get_option('output_filename')
        # output_dir = self.get_parser_settings().get('output_dir', None)
        output_dir = self.node.get_option('output_dir')
        #output_dir = self.node.get_option('output_dir')
        #output_dir = self.node.get_incoming(node_class=str, link_label_filter='output_dir').one().node

        

        # Check that folder content is as expected
        # This simple check makes sure that the expected output file 
        # file.log is among the files retrieved from the computer where 
        # the calculation was run. 
        files_retrieved = self.retrieved.list_object_names()
        files_expected = [] #[output_filename]
        if 'output_files' in self.node.inputs:
            for name in self.node.inputs.output_files:
                files_expected.extend([str(name)])
        # Note: set(A) <= set(B) checks whether A is a subset of B
        # if not set(files_expected) <= set(files_retrieved):
        #     self.logger.error(f"Found files '{files_retrieved}', expected to find '{files_expected}'")
        #     # AiiDA stores the exit code returned by the parse method on the 
        #     # calculation node that is being parsed.
        #     # exit code is defined in CalcJob.
        #     return self.exit_codes.ERROR_MISSING_OUTPUT_FILES 

        for file in files_expected:
            if file not in files_retrieved:
                self.logger.error(f"User defined output file '{file}' not in "
                                f"list of retrieved files '{files_retrieved}'")
                return self.exit_codes.ERROR_UNTRACKED_OUTPUT_FILES

        # add output file
        # simply passing along the output file as a SinglefileData node.
        # instead of just returning the file you may want to parse it 
        # e.g. to a python dictionary (Dict node) to make the results easily 
        # searchable.

        self.logger.info(f"Parsing '{output_filename}'")
        with self.retrieved.open(output_filename, 'rb') as handle:
            output_node = SinglefileData(file=handle)

        # the out() method is used return the output file as the exec output 
        # of the calculation: The first argument is the name to be used as 
        # the label for the link that connects the calculation and data node. 
        # The second argument is the node that should be recorded as an output.
        
        self.out('log', output_node)

        for index, thing in enumerate(files_expected):
            self.logger.info(f"Parsing '{thing}'")
            with self.retrieved.open(thing, 'rb') as handle:
                output_node = SinglefileData(file=handle, filename=thing)
            self.out(self.format_link_label(thing), output_node)

        for index, thing in enumerate(files_retrieved):
            self.logger.info(f"Parsing '{thing}'")
            file_path = os.path.join(output_dir, thing)
            try:
                with self.retrieved.open(thing, 'rb') as handle:
                    with open(file_path, 'wb') as f_out:
                        while True:
                            # Read a chunk of the input file
                            chunk = handle.read(1024)
                            if not chunk:
                                # End of file reached
                                break
                            # Write the chunk to the output file
                            f_out.write(chunk)
            except UnicodeDecodeError:
                with self.retrieved.open(thing, 'r') as handle:
                    f = open(file_path, 'w')
                    for line in handle.read():
                        f.write(line)
                    f.close()

        # Note: The outputs and their types need to match those from the 
        # process specification of the corresponding CalcJob 
        # (or an exception will be raised).

        # In order to request automatic parsing of a CalcJob 
        # (once it has finished), users can set the 
        # metadata.options.parser_name input when launching the job. 
        # If a particular parser should be used by default, the CalcJob define 
        # method can set a default value for the parser name.
        # `@classmethod
        #   def define(cls, spec):
        #       ...
        #       spec.inputs['metadata']['options']['parser_name'].default = \
        #           'general-MD'`
        #
        # Note that the default is not set to the Parser class itself, 
        # but to the entry point string under which the parser class is 
        # registered. 


        return ExitCode(0)
    

    @staticmethod
    def format_link_label(filename: str) -> str:
        """
        Modified from: https://github.com/sphuber/aiida-shell/blob/master/src/aiida_shell/parsers/shell.py
        Format the link label from a given filename.
        Valid link labels can only contain alphanumeric characters and
            underscores, without consecutive underscores. So all characters 
            that are not alphanumeric or an underscore are converted to 
            underscores, where consecutive underscores are merged into one.
        Additional: Label cannot start with a number or underscore.
        :param filename: The filename.
        :returns: The link label.
        """
        filename = filename.split('/')[-1]
        if filename[0].isdigit() or filename[0] == '_':
            filename = 'output_files_' + filename
        alphanumeric = re.sub('[^0-9a-zA-Z_]+', '_', filename)
        link_label = re.sub('_[_]+', '_', alphanumeric)

        return link_label
