import os

class GenOrca:
    # Class for storing information for Gaussian jobs and generating scripts
    def __init__(self, 
                 directory, 
                 name,
                 xyz=None,
                 opt=False,
                 freeze_atoms=[],
                 freeze_bonds=[],
                 ts=False,
                 freq=False,
                 irc=False,
                 irc_direction=None,
                 irc_opt=False,
                 charge=None,
                 multiplicity=1,
                 nprocs=8, 
                 memory=8,
                 functional="B3LYP",
                 basis=None,
                 dispersion_method=None,
                 solvent=None,
                 route=None,
                 output_name=None,
                 geom_chk=False,
                 old_checkpoint=None,
                 no_freeze=False,
                 nbo=False
                 ) -> None:
        self.directory = directory
        self.name = name
        self.geom_chk = geom_chk
        if not xyz:
            self.xyz = f"{name}.xyz"
        else:
            self.xyz = xyz
        self.opt = opt
        self.freeze_atoms=freeze_atoms
        self.freeze_bonds = freeze_bonds
        self.ts = ts
        self.freq = freq
        self.irc = irc
        self.irc_direction = irc_direction
        self.irc_opt = irc_opt
        self.charge = charge
        self.multiplicity = multiplicity
        self.nprocs = nprocs
        self.memory = memory
        self.functional = functional
        self.basis = basis
        self.dispersion_method = dispersion_method
        self.solvent = solvent
        self.old_checkpoint = old_checkpoint
        self.route = route
        self.no_freeze = no_freeze
        self.nbo = nbo
        if not self.route:
            self.route = []
        self.title = self.gen_title()
        if not output_name:
            self.output_name = self.gen_output_name()
        else:
            self.output_name = output_name
        if not self.charge:
            self.get_charge()
        self.read_xyz()

    def gen_output(self):
        full_path = os.path.join(self.directory,self.output_name)
        with open(f"{full_path}.inp", "w") as f:
            f.write(f"{self.gen_input_line()}\n")
            f.write(f"%maxcore {self.memory * 1000}\n")
            f.write("%pal\n")
            f.write(f"nprocs {self.nprocs}\n")
            f.write("end\n")
            f.write(f"# {self.title}\n")
            f.write("\n")
            f.write(f"*xyz {self.charge} {self.multiplicity}\n")
            for line in self.xyz_lines:
                f.write(line)
            f.write("\n")
            f.write("*")
            f.close()

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
            f.write("\n")
            f.write(f"name={self.output_name}\n")
            f.write("\n")
            f.write("# Initialise the environment modules\n")
            f.write(". /etc/profile.d/modules.sh\n")
            f.write("\n")
            f.write("# Export environment variables\n")
            f.write('export PATH="/exports/eddie3_homes_local/dwillco2/openmpi502/bin:$PATH"\n')
            f.write('export LD_LIBRARY_PATH="/exports/eddie3_homes_local/dwillco2/openmpi502/lib:$LD_LIBRARY_PATH"\n')
            f.write('ORCA_DIR="/exports/csce/eddie/chem/groups/Thomas"\n')
            f.write('export PATH="$ORCA_DIR/orca:$PATH"\n')
            f.write('export LD_LIBRARY_PATH="$ORCA_DIR/orca:$LD_LIBRARY_PATH"\n')
            f.write('SCR=${TMPDIR}\n')
            f.write('ORIG=`pwd`\n')
            f.write('cp ${name}.inp $SCR\n')
            f.write('cd $SCR\n')
            f.write("\n")
            f.write("# Run the program\n")
            f.write('$ORCA_DIR/orca/orca "${name}.inp" > "${ORIG}/${name}.out"\n')
            f.write('cp * $ORIG\n')
            f.write('cd $ORIG\n')
            f.close()
    
    def gen_output_name(self):
        if self.irc:
            if self.solvent:
                return f"irc_{self.irc_direction}_{self.name}_{self.solvent}"
            else:
                return f"irc_{self.irc_direction}_{self.name}"
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
                return f"{self.functional.lower()}_{self.name}_{self.solvent}"
            else:
                return f"{self.functional.lower()}_{self.name}"

    def gen_input_line(self):
        ri_input = False
        input_line = "! "
        if "RI-" in self.functional:
            ri_input = True
        input_line += f"{self.functional} "
        if self.dispersion_method:
            input_line += f"{self.dispersion_method} "
        if self.basis:
            input_line += f"{self.basis} "
        if ri_input:
            input_line += "RIJCOSX AutoAux DefGrid2 "
        input_line += "TightSCF "
        if self.opt:
            input_line += "Opt "
        if self.freq:
            input_line += "Freq "
        return input_line
    
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
            
    def gen_title(self):
        if self.irc:
            return f"{self.name} irc {self.irc_direction}"
        elif self.opt:
            return f"{self.name} opt"
        elif self.freq:
            return f"{self.name} freq"
        else:
            return f"{self.name}"
    
    def read_xyz(self):
        xyz_directory = os.path.join(self.directory,self.xyz)
        with open(xyz_directory, "r") as f:
            xyz_lines = f.readlines()
            self.xyz_lines = xyz_lines[2:]
            f.close()
 