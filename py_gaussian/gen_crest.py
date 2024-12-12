import os
from .gen_job import GenJob

CREST_SETUP = r"""flight env activate gridware
module load apps/crest
export TMPDIR=/tmp/users/drw4001/$SLURM_JOB_ID
export SCR=${TMPDIR}
mkdir -p ${TMPDIR}
# Copy the input data to the directory pointed to by the $TMPDIR variable
cp ${SLURM_SUBMIT_DIR}/* ${TMPDIR}
# Go to the $TMPDIR directory
cd $TMPDIR
"""

class GenCrest (GenJob):
        # Class for storing information for CREST jobs and generating scripts
    def __init__(self, 
                 directory, 
                 name,
                 xyz=None,
                 freeze_atoms=[],
                 charge=None,
                 multiplicity=1,
                 nprocs=8, 
                 memory=8000,
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
                 solvent)
        self.freeze_atoms = freeze_atoms
        if not output_name:
            self.output_name = self.gen_output_name()
        else:
            self.output_name = output_name

    def gen_bash(self, runtime=48):
        template = self.gen_bash_template(runtime=runtime)
        full_path = os.path.join(self.directory,self.output_name)
        with open(f"{full_path}.sh", "w") as f:
            f.write(template)
            f.write(CREST_SETUP)
            f.write(f"INPUT={self.xyz}\n")
            if self.freeze_atoms:
                f.write(f"crest $INPUT --constrain {','.join(map(str, self.freeze_atoms))}\n")
                if self.solvent:
                    f.write(f"crest $INPUT --gfn2 --cinp .xcontrol.sample --chrg {self.charge} --alpb {self.solvent} -T {self.nprocs} > $SLURM_SUBMIT_DIR/{self.output_name}.out\n")
                else:
                    f.write(f"crest $INPUT --gfn2 --cinp .xcontrol.sample --chrg {self.charge} -T {self.nprocs} > $SLURM_SUBMIT_DIR/{self.output_name}.out\n")
            elif self.solvent:
                f.write(f"crest $INPUT --gfn2 --chrg {self.charge} --alpb {self.solvent} -T {self.nprocs} > $SLURM_SUBMIT_DIR/{self.output_name}.out\n")
            else:
                f.write(f"crest $INPUT --gfn2 --chrg {self.charge} -T {self.nprocs} > $SLURM_SUBMIT_DIR/{self.output_name}.out\n")
            f.write(f"cp -r $TMPDIR/* $SLURM_SUBMIT_DIR/\n")
            f.write(f"rm -rf $TMPDIR\n")
            f.close()

    def gen_output_name(self):
        if self.solvent:
            return f"crest_{self.name}_{self.solvent}"
        else:
            return f"crest_{self.name}"

    
    