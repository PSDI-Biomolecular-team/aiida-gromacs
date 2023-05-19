"""
aiida_gromacs

A workflow for setting up basic molecular dynamics simulations.
"""
from aiida.engine import ToContext, WorkChain
from aiida.orm import Code, SinglefileData
from aiida.plugins.factories import CalculationFactory, DataFactory

from aiida_gromacs import helpers

Pdb2gmxCalculation = CalculationFactory("gromacs.pdb2gmx")
EditconfCalculation = CalculationFactory("gromacs.editconf")
SolvateCalculation = CalculationFactory("gromacs.solvate")
GromppCalculation = CalculationFactory("gromacs.grompp")
GenionCalculation = CalculationFactory("gromacs.genion")
MdrunCalculation = CalculationFactory("gromacs.mdrun")

Pdb2gmxParameters = DataFactory("gromacs.pdb2gmx")
EditconfParameters = DataFactory("gromacs.editconf")
SolvateParameters = DataFactory("gromacs.solvate")
GromppParameters = DataFactory("gromacs.grompp")
GenionParameters = DataFactory("gromacs.genion")
MdrunParameters = DataFactory("gromacs.mdrun")


class SetupWorkChain(WorkChain):
    """WorkChain for setting up a gromacs simulation automatically."""

    @classmethod
    def define(cls, spec):
        """Specify workflow recipe."""
        super().define(spec)
        spec.input("local_code", valid_type=Code)
        spec.input("remote_code", required=False, valid_type=Code)
        spec.input("pdbfile", valid_type=SinglefileData, help="Input structure.")
        spec.input(
            "ionsmdp", valid_type=SinglefileData, help="MD parameters for adding ions."
        )
        spec.input(
            "minmdp", valid_type=SinglefileData, help="MD parameters for minimisation."
        )
        spec.input(
            "nvtmdp",
            valid_type=SinglefileData,
            help="MD parameters for NVT equilibration.",
        )
        spec.input(
            "nptmdp",
            valid_type=SinglefileData,
            help="MD parameters for NPT equilibration.",
        )
        spec.input(
            "prodmdp",
            valid_type=SinglefileData,
            help="MD parameters for production run.",
        )
        spec.input(
            "pdb2gmxparameters",
            valid_type=Pdb2gmxParameters,
            help="Command line parameters for gmx pdb2gmx",
        )
        spec.input(
            "editconfparameters",
            valid_type=EditconfParameters,
            help="Command line parameters for gmx editconf",
        )
        spec.input(
            "solvateparameters",
            valid_type=SolvateParameters,
            help="Command line parameters for gmx solvate",
        )
        spec.input(
            "gromppionsparameters",
            valid_type=GromppParameters,
            help="Command line parameters for gmx grompp",
        )
        spec.input(
            "genionparameters",
            valid_type=GenionParameters,
            help="Command line parameters for gmx genion",
        )
        spec.input(
            "gromppminparameters",
            valid_type=GromppParameters,
            help="Command line parameters for gmx grompp minimisation run",
        )
        spec.input(
            "minimiseparameters",
            valid_type=MdrunParameters,
            help="Command line parameters for gmx mdrun minimisation run",
        )
        spec.input(
            "gromppnvtparameters",
            valid_type=GromppParameters,
            help="Command line parameters for gmx grompp nvt equilibration run",
        )
        spec.input(
            "nvtparameters",
            valid_type=MdrunParameters,
            help="Command line parameters for gmx mdrun nvt equilibration run",
        )
        spec.input(
            "gromppnptparameters",
            valid_type=GromppParameters,
            help="Command line parameters for gmx grompp npt equilibration run",
        )
        spec.input(
            "nptparameters",
            valid_type=MdrunParameters,
            help="Command line parameters for gmx mdrun npt equilibration run",
        )
        spec.input(
            "gromppprodparameters",
            valid_type=GromppParameters,
            help="Command line parameters for gmx grompp production run",
        )
        spec.input(
            "mdrunparameters",
            valid_type=MdrunParameters,
            help="Command line parameters for gmx mdrun production run",
        )

        spec.outline(
            cls.pdb2gmx,
            cls.editconf,
            cls.solvate,
            cls.gromppions,
            cls.genion,
            cls.gromppmin,
            cls.minimise,
            cls.gromppnvt,
            cls.nvtequilibrate,
            cls.gromppnpt,
            cls.nptequilibrate,
            cls.gromppprod,
            cls.prodmd,
            cls.result,
        )

        spec.output("result")

    def pdb2gmx(self):
        """Convert PDB file to forcefield compliant GRO file"""
        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.pdb2gmxparameters,
            "pdbfile": self.inputs.pdbfile,
            "metadata": {
                "description": "convert pdb file to gromacs gro format.",
            },
        }

        future = self.submit(Pdb2gmxCalculation, **inputs)

        return ToContext(pdb2gmx=future)

    def editconf(self):
        """Add simulation box to GRO file."""

        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.editconfparameters,
            "grofile": self.ctx.pdb2gmx.outputs.grofile,
            "metadata": {
                "description": "add simulation box parameters to gro file.",
            },
        }

        future = self.submit(EditconfCalculation, **inputs)

        return ToContext(editconf=future)

    def solvate(self):
        """Add solvent to GRO file."""
        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.solvateparameters,
            "grofile": self.ctx.editconf.outputs.grofile,
            "topfile": self.ctx.pdb2gmx.outputs.topfile,
            "metadata": {
                "description": "add solvent to simulation box.",
            },
        }

        future = self.submit(SolvateCalculation, **inputs)

        return ToContext(solvate=future)

    def gromppions(self):
        """Create a tpr for adding ions."""
        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.gromppionsparameters,
            "mdpfile": self.inputs.ionsmdp,
            "grofile": self.ctx.solvate.outputs.grofile,
            "topfile": self.ctx.solvate.outputs.topfile,
            "metadata": {
                "description": "prepare the tpr for adding ions.",
            },
        }

        future = self.submit(GromppCalculation, **inputs)

        return ToContext(gromppions=future)

    def genion(self):
        """Add ions to system to balance charge."""

        # Sort this out.
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point="bash", computer=computer)
        inputs = {
            "code": gromacs_code,
            "parameters": self.inputs.genionparameters,
            "tprfile": self.ctx.gromppions.outputs.tprfile,
            "topfile": self.ctx.solvate.outputs.topfile,
            "metadata": {
                "description": "add ions to simulation box.",
            },
        }

        future = self.submit(GenionCalculation, **inputs)

        return ToContext(genion=future)

    def gromppmin(self):
        """Create a tpr for minimisation."""
        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.gromppminparameters,
            "mdpfile": self.inputs.minmdp,
            "grofile": self.ctx.genion.outputs.grofile,
            "topfile": self.ctx.genion.outputs.topfile,
            "metadata": {
                "description": "prepare the tpr for minimisation run.",
            },
        }

        future = self.submit(GromppCalculation, **inputs)

        return ToContext(gromppmin=future)

    def minimise(self):
        """Minimise system."""
        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.minimiseparameters,
            "tprfile": self.ctx.gromppmin.outputs.tprfile,
            "metadata": {
                "description": "minimise system.",
            },
        }

        future = self.submit(MdrunCalculation, **inputs)

        return ToContext(minimise=future)

    def gromppnvt(self):
        """Create a tpr for NVT equilibration."""
        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.gromppnvtparameters,
            "mdpfile": self.inputs.nvtmdp,
            "grofile": self.ctx.minimise.outputs.grofile,
            "topfile": self.ctx.genion.outputs.topfile,
            "itpfile": self.ctx.pdb2gmx.outputs.itpfile,
            "metadata": {
                "description": "prepare the tpr for NVT equlibration.",
            },
        }

        future = self.submit(GromppCalculation, **inputs)

        return ToContext(gromppnvt=future)

    def nvtequilibrate(self):
        """NVT Equilibration of system."""
        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.nvtparameters,
            "tprfile": self.ctx.gromppnvt.outputs.tprfile,
            "metadata": {
                "description": "NVT equilibrate system.",
            },
        }

        future = self.submit(MdrunCalculation, **inputs)

        return ToContext(nvtequilibrate=future)

    def gromppnpt(self):
        """Create a tpr for NPT equilibration."""
        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.gromppnptparameters,
            "mdpfile": self.inputs.nptmdp,
            "grofile": self.ctx.nvtequilibrate.outputs.grofile,
            "topfile": self.ctx.genion.outputs.topfile,
            "itpfile": self.ctx.pdb2gmx.outputs.itpfile,
            "metadata": {
                "description": "prepare the tpr for NPT equlibration.",
            },
        }

        future = self.submit(GromppCalculation, **inputs)

        return ToContext(gromppnpt=future)

    def nptequilibrate(self):
        """NPT Equilibration of system system."""
        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.nptparameters,
            "tprfile": self.ctx.gromppnpt.outputs.tprfile,
            "metadata": {
                "description": "NPT equilibrate system.",
            },
        }

        future = self.submit(MdrunCalculation, **inputs)

        return ToContext(nptequilibrate=future)

    def gromppprod(self):
        """Create a tpr for production run."""
        inputs = {
            "code": self.inputs.local_code,
            "parameters": self.inputs.gromppprodparameters,
            "mdpfile": self.inputs.prodmdp,
            "grofile": self.ctx.nptequilibrate.outputs.grofile,
            "topfile": self.ctx.genion.outputs.topfile,
            "itpfile": self.ctx.pdb2gmx.outputs.itpfile,
            "metadata": {
                "description": "prepare the tpr for production run.",
            },
        }

        future = self.submit(GromppCalculation, **inputs)

        return ToContext(gromppprod=future)

    def prodmd(self):
        """Run production MD"""

        if "remote_code" in self.inputs:
            code = self.inputs.remote_code

        else:
            code = self.inputs.local_code

        inputs = {
            "code": code,
            "parameters": self.inputs.mdrunparameters,
            "tprfile": self.ctx.gromppprod.outputs.tprfile,
            "metadata": {
                "description": "Production MD.",
            },
        }

        future = self.submit(MdrunCalculation, **inputs)

        return ToContext(prodmd=future)

    def result(self):
        """Results"""
        self.out("result", self.ctx.prodmd.outputs.trrfile)
