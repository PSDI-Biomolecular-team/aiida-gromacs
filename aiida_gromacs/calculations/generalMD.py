
from aiida.common import datastructures
from aiida.engine import CalcJob
from aiida.orm import SinglefileData, Str, List

class GeneralCalculation(CalcJob):
    """
    AiiDA calculation plugin wrapping an executable with user defined
    input and output files.
    """
    @classmethod
    def define(cls, spec):
        """
        Define inputs and outputs of the calculation.
            The define method tells AiiDA which inputs the CalcJob 
        expects and which outputs it produces.
            This is done through an instance of the CalcJobProcessSpec class,
        which is passed as the spec argument to the define method.
        """
        # yapf: disable

        # line below calls the define method of the CalcJob parent class. 
        # This necessary step defines the inputs and outputs that are common 
        # to all CalcJobs.
        super(GeneralCalculation, cls).define(spec)

        # new ports
        # use the input() method to define our two input files
        # When using SinglefileData, AiiDA keeps track of the inputs as files.
        spec.input('command',
                valid_type=Str, required=False,
                help='The command used to execute the job.')
        
        # spec.input('output_dir', valid_type=str, 
        #         help='The directory where output files will be saved.')
        
        # dynamic dict of inputs instead of specifying each one.
        spec.input_namespace(
            'input_files',
            valid_type=SinglefileData,
            required=False,
            help='Dictionary of input files.',
            dynamic=True, # can take number of values unknown at time of definition
        )

        # use output() to define the only output of the calculation 
        # with the label execute.
        # AiiDA will attach the outputs defined here to a (successfully) 
        # finished calculation using the link label provided.
        spec.output('log', valid_type=SinglefileData, required=False,
                help='link to the default file.log.')
        
        # set the list of output file names as an input so that it can be
        # iterated over in the parser later.
        spec.input('output_files', valid_type=List, required=False, 
                   help='List of output file names.')

        # define the schema for metadata.options
        spec.input('metadata.options.output_filename', valid_type=str, 
                default='file.log', help='name of file produced by default.')
        spec.input('metadata.options.output_dir', valid_type=str, 
                help='The directory where output files will be saved '
                    'when parsed.')
        
        # spec.input('metadata.options.parser_name', 
        #         valid_type=Str, required=False, default=Str('general-MD'),
        #         help='The name of the parser to use.')
        
        # set a few default options below. These options have already 
        # been defined on the spec by the super().define(spec) call, 
        # and they can be accessed through the inputs attribute, 
        # which behaves like a dictionary.
        # set name of output file.
        
        # set the computational resources used for this calculation.
        spec.inputs['metadata']['options']['resources'].default = {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 1,
                }
        # set name of the parser.
        spec.inputs['metadata']['options']['parser_name'].default = \
                'general-MD'
        
        # ensure code is set
        spec.inputs['code'].required = True

        
        # IMPORTANT: 
        # Use spec.outputs.dynamic = True to make the entire output namespace 
        # fully dynamic. This means any number of output files 
        # can be linked to a node.
        spec.outputs.dynamic = True
        spec.inputs.dynamic = True
        spec.inputs['metadata']['options'].dynamic = True

        # Used by parsers, which communicate errors through exit codes.
        # An exit_code defines:
        #   -   an exit status (a positive integer, following the Exit code 
        #       conventions),
        #   -   a label that can be used to reference the code in the parse 
        #       method self.exit_codes property, and
        #   -   a message that provides a more detailed description of the 
        #       problem.
        spec.exit_code(300, 'ERROR_MISSING_OUTPUT_FILES',
                message='Calculation did not produce all expected '\
                        'output files.')
        spec.exit_code(301, 'ERROR_UNTRACKED_OUTPUT_FILES',
                message='Specified output file not produced by command.')

        # There is no return statement in define: the define method 
        # directly modifies the spec object it receives.

        # One more input required by any CalcJob is which external 
        # executable to use.
        # External executables are represented by AbstractCode instances 
        # that contain information about the computer they reside on, 
        # their path in the file system and more. 
        # They are passed to a CalcJob via the code input, 
        # which is defined in the CalcJob base class, so you don’t have to.
        # spec.input('code', valid_type=orm.AbstractCode, 
        #     help='The `Code` to use for this job.')


    def prepare_for_submission(self, folder):
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin 
            should temporarily place all files needed by the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance

        The prepare_for_submission() method has two jobs: 
            Creating the input files in the format the external 
            code expects and returning a CalcInfo object that contains 
            instructions for the AiiDA engine on how the code should be run.

        All inputs provided to the calculation are validated against 
            the spec before prepare_for_submission() is called. 
            Therefore, when accessing the inputs attribute, you can 
            safely assume that all required inputs have been set and 
            that all inputs have a valid type.
        """

        # start by creating a CodeInfo object that lets AiiDA know how 
        # to run the code,
        codeinfo = datastructures.CodeInfo()


        
        # split strings in command
        codeinfo.cmdline_params = str(self.inputs.command.value).split()
        if self.inputs.code.label == 'bash':
            codeinfo.cmdline_params = ['-c', self.inputs.command.value]

        # and the UUID of the AbstractCode to run:
        codeinfo.code_uuid = self.inputs.code.uuid

        # Since anything that writes directly to standard output, we redirect 
        # standard output to the specified output filename.
        codeinfo.stdout_name = self.metadata.options.output_filename
        # output_dir = self.metadata.options.output_dir

        #calcinfo.parser_args = {'output_dir': self.inputs.output_dir}

        # Prepare a `CalcInfo` to be returned to the engine
        # create a CalcInfo object that lets AiiDA know 
        # which files to copy back and forth. 
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]

        # get a list of input files to be copied to remote.
        copy_list = [] # used below
        if 'input_files' in self.inputs:
            for name, obj in self.inputs.input_files.items():
                copy_list.append((obj.uuid, obj.filename, obj.filename))

        # The input files are already stored 
        # in the AiiDA file repository and we can use the local_copy_list 
        # to pass them along.
        calcinfo.local_copy_list = copy_list

        # In other use cases, may need to create new files on the fly. 
        # So use the folder argument of prepare_for_submission():
        # `with folder.open("filename", 'w') as handle:
        #   handle.write("file content")`
        # Any files and directories created in this sandbox folder will 
        # automatically be transferred to the compute resource where the 
        # actual calculation takes place.

        # The retrieve_list tells the engine which files to retrieve 
        # from the directory where the job ran after it has finished. 
        # All files listed here will be stored in a FolderData node 
        # that is attached as an output node to the calculation with 
        # the label retrieved.
        retrieve_list = [self.metadata.options.output_filename]
        if 'output_files' in self.inputs: # check there are output files.
            for name in self.inputs.output_files:
                retrieve_list.append(str(name)) # save output filename to list
        calcinfo.retrieve_list = retrieve_list
        calcinfo.retrieve_temporary_list = retrieve_list

        return calcinfo
