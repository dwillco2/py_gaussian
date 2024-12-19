import os
from .gen_job import GenJob

ORCA_BASH = r"""flight env activate gridware
module load apps/orca
module load mpi/openmpi
export TMPDIR=/tmp/users/drw4001/$SLURM_JOB_ID
export SCR=${TMPDIR}
mkdir -p ${TMPDIR}
cp ${SLURM_SUBMIT_DIR}/* ${TMPDIR}
cd $TMPDIR
"""

class GenOrca(GenJob):
    # Class for storing information for Gaussian jobs and generating scripts
    def __init__(self, 
                 directory, 
                 name,
                 xyz=None,
                 opt=False,
                 goat=False,
                 freeze_atoms=[],
                 freeze_bonds=[],
                 ts=False,
                 freq=False,
                 scan=False,
                 scan_list=None, # [atom_0, atom_1, d_0, d_1, nsteps] 
                 irc=False,
                 irc_opt=False,
                 charge=None,
                 multiplicity=1,
                 nprocs=8, 
                 memory=8000,
                 functional="R2SCAN-3c",
                 basis=None,
                 dispersion_method=None,
                 solvent=None,
                 output_name=None,
                 ) -> None:
        super().__init__(directory, 
                 name,
                 xyz,
                 charge,
                 multiplicity,
                 nprocs, 
                 memory,
                 solvent,)
        self.opt = opt
        self.freeze_atoms=freeze_atoms
        self.freeze_bonds = freeze_bonds
        self.ts = ts
        self.freq = freq
        self.scan = scan
        self.scan_list = scan_list
        self.irc = irc
        self.irc_opt = irc_opt
        self.goat = goat
        if self.goat:
            self.functional = "XTB"
        else:
            self.functional = functional
        self.basis = basis
        self.dispersion_method = dispersion_method
        self.title = self.gen_title()
        if not output_name:
            self.output_name = self.gen_output_name()
        else:
            self.output_name = output_name

    def gen_output(self):
        full_path = os.path.join(self.directory,self.output_name)
        with open(f"{full_path}.inp", "w") as f:
            f.write(f"{self.gen_input_line()}\n")
            f.write(f"%maxcore {self.memory}\n")
            f.write("%pal\n")
            f.write(f"nprocs {self.nprocs}\n")
            f.write("end\n")
            if self.scan:
                f.write("%geom Scan\n")
                f.write(f"B {self.scan_list[0]} {self.scan_list[1]} = {self.scan_list[2]}, {self.scan_list[3]}, {self.scan_list[4]}\n")
                f.write("end\n")
                f.write("end\n")
            if self.ts:
                f.write("%geom\n")
                f.write("Calc_Hess true\n")
                f.write("NumHess true\n")
                f.write("Recalc_Hess 5\n")
                f.write("end\n")
            if self.irc:
                f.write("%irc\n")
                f.write(" MAXITER 30\n")
                f.write("end\n") 
            f.write(f"# {self.title}\n")
            f.write("\n")
            f.write(f"*XYZFILE {self.charge} {self.multiplicity} {self.xyz}\n")
            f.close()

    def gen_bash(self, runtime=48):
        template = self.gen_bash_template(runtime=runtime)
        full_path = os.path.join(self.directory,self.output_name)
        with open(f"{full_path}.sh", "w") as f:
            f.write(template)
            f.write(ORCA_BASH)
            f.write(f"$ORCADIR/orca {self.output_name}.inp > $SLURM_SUBMIT_DIR/{self.output_name}.out\n")
            f.write("cp -r $TMPDIR/* $SLURM_SUBMIT_DIR/\n")
            f.write("rm -rf $TMPDIR\n")
            f.close()

    
    def gen_output_name(self):
        if self.goat:
            if self.solvent:
                return f"goat_{self.name}_{self.solvent}"
            else:
                return f"goat_{self.name}"
        elif self.scan:
            if self.solvent:
                return f"scan_{self.name}_{self.solvent}"
            else:
                return f"scan_{self.name}"
        elif self.irc:
            if self.solvent:
                return f"irc_{self.name}_{self.solvent}"
            else:
                return f"irc_{self.name}"
        elif self.opt and self.ts:
            if self.solvent:
                return f"opt_ts_{self.name}_{self.solvent}"
            else:
                return f"opt_ts_{self.name}"
        elif self.opt and not self.ts:
            if self.solvent:
                return f"opt_{self.name}_{self.solvent}"
            else:
                return f"opt_{self.name}"
        elif self.freq:
            if self.solvent:
                return f"freq_{self.name}_{self.solvent}"
            else:
                return f"freq_{self.name}"
        else:
            if self.solvent:
                return f"sp_{self.name}_{self.solvent}"
            else:
                return f"sp_{self.name}"

    def gen_input_line(self):
        ri_input = False
        input_line = "! "
        if self.goat:
            input_line += "GOAT "
        if "RI-" in self.functional:
            ri_input = True
        if self.irc:
            input_line += "IRC "
        input_line += f"{self.functional} "
        if self.dispersion_method:
            input_line += f"{self.dispersion_method} "
        if self.basis:
            input_line += f"{self.basis} "
        if ri_input:
            input_line += "RIJCOSX AutoAux DefGrid2 "
        if self.opt or self.scan:
            if not self.ts:
                input_line += "Opt "
        if self.opt and self.ts:
            input_line += "OptTS "
        if self.freq and not self.ts:
            input_line += "Freq "
        if self.freq and self.ts:
            input_line += "NumFreq "
        if not self.goat:
            input_line += "TightSCF "
        return input_line
    
            
    def gen_title(self):
        if self.irc:
            return f"{self.name} irc"
        elif self.opt:
            return f"{self.name} opt"
        elif self.freq:
            return f"{self.name} freq"
        else:
            return f"{self.name}"
    
 