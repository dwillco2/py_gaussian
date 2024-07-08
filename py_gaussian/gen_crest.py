import os

class GenCrest:
        # Class for storing information for CREST jobs and generating scripts
    def __init__(self, 
                 directory, 
                 name,
                 xyz=None,
                 opt=False,
                 freeze_atoms=[],
                 charge=None,
                 multiplicity=1,
                 nprocs=8, 
                 memory=8,
                 solvent=None,
                 output_name=None,
                 ) -> None:
        self.directory = directory
        self.name = name
        self.xyz = xyz
        self.opt = opt
        self.freeze_atoms=freeze_atoms
        self.charge = charge
        self.multiplicity = multiplicity
        self.nprocs = nprocs
        self.memory = memory
        self.solvent = solvent
        if not output_name:
            self.output_name = self.gen_output_name()
        else:
            self.output_name = output_name
        if not self.charge:
            self.get_charge()

    def gen_bash(self, runtime=48):
        full_path = os.path.join(self.directory,self.output_name)
        with open(f"{full_path}.bash", "w") as f:
            f.write("#!/bin/sh\n")
            f.write("# Grid Engine options (lines prefixed with #$)\n")
            f.write(f"#$ -N {self.output_name}\n")
            f.write("#$ -cwd\n")
            f.write(f"#$ -l h_rt={runtime}:00:00\n")
            f.write(f"#$ -l h_vmem={self.memory}G\n")
            f.write(f"#$ -pe sharedmem {self.nprocs}\n")
            f.write(f"#$ -l rl9=true\n")
            f.write("\n")
            f.write("# Initialise the environment modules\n")
            f.write(". /etc/profile.d/modules.sh\n")
            f.write("\n")
            f.write("# Export environment variables\n")
            f.write("module load anaconda\n")
            f.write("source activate crest2\n")
            f.write(f"INPUT={self.xyz}\n")
            f.write("\n")
            f.write("# Run the program\n")
            if self.solvent:
                f.write(f"crest $INPUT --gfn2 --chrg {self.charge} --alpb {self.solvent} -T {self.nprocs}\n")
            else:
                f.write(f"crest $INPUT --gfn2 --chrg {self.charge} -T {self.nprocs}\n")
            f.close()

    def gen_output_name(self):
        if self.solvent:
            return f"crest_{self.name}_{self.solvent}"
        else:
            return f"crest_{self.name}"

    def get_charge(self):
        for filename in os.listdir(self.directory):
            i = os.path.join(self.directory, filename)
            if os.path.isfile(i) and ".CHRG" in i:
                with open(i, "r") as f:
                    lines = f.readlines()
                    self.charge = int(lines[0])
                    f.close()
                break
        else:
            self.charge = 0
            print("charge file not found, charge = 0")
    
    