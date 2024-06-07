import os

class GenGaussian:
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
                 basis="Def2SVP",
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
        if not xyz and not self.geom_chk:
            self.xyz = f"{name}.xyz"
        elif self.geom_chk:
            self.xyz = None
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
        with open(f"{full_path}.gjf", "w") as f:
            if self.irc:
                if self.old_checkpoint:
                    old_chk = self.old_checkpoint
                else:
                    raise Exception("old checkpoint file required for irc calculation")
                if ".chk" in old_chk:
                    f.write(f"%oldchk={old_chk} \n")
                else:
                    f.write(f"%oldchk={old_chk}.chk \n")
            elif self.freq and not self.opt and self.geom_chk:
                if self.old_checkpoint:
                    old_chk = self.old_checkpoint
                elif self.ts:
                    old_chk = self.output_name.replace("freq","opt_ts")
                else:
                    old_chk = self.output_name.replace("freq","opt")
                if ".chk" in old_chk:
                    f.write(f"%oldchk={old_chk} \n")
                else:
                    f.write(f"%oldchk={old_chk}.chk \n")
            elif self.opt and self.geom_chk:
                if self.old_checkpoint:
                    old_chk = self.old_checkpoint
                else:
                    raise Exception("geom=allcheck given, but no checkpoint file provided")
                if ".chk" in old_chk:
                    f.write(f"%oldchk={old_chk} \n")
                else:
                    f.write(f"%oldchk={old_chk}.chk \n")
            f.write(f"%chk={self.output_name}.chk \n")
            f.write(f"%nproc={self.nprocs} \n")
            f.write(f"%mem={self.memory}GB \n")
            f.write(f"{self.gen_input_line()}\n")
            f.write("\n")
            f.write(f"{self.title}\n")
            f.write("\n")
            if not self.geom_chk:
                f.write(f"{self.charge} {self.multiplicity}\n")
                for line in self.xyz_lines:
                    f.write(line)
                f.write("\n")
            if len(self.freeze_atoms) > 0:
                for atom in self.freeze_atoms:
                    f.write(f"{atom} F\n")
            if len(self.freeze_bonds) > 0:
                for bond in self.freeze_bonds:
                    f.write(f"{bond[0]} {bond[1]} F\n")
            if self.nbo:
                f.write("$nbo bndidx $end\n")
            f.write("\n")
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
            f.write(f"input={self.output_name}.gjf\n")
            f.write("\n")
            f.write("# Initialise the environment modules\n")
            f.write(". /etc/profile.d/modules.sh\n")
            f.write("\n")
            f.write("# Export environment variables\n")
            f.write("export g16root=/exports/applications/apps/community/chem\n")
            f.write("export GAUSS_SCRDIR=$TMPDIR\n")
            f.write("source $g16root/g16/bsd/g16.profile\n")
            f.write("\n")
            f.write("# Run the program\n")
            f.write("$g16root/g16/g16 $input\n")
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

    def gen_input_line(self):
        input_line = "#p "
        mod_redun = False
        if len(self.freeze_atoms) > 0 or len(self.freeze_bonds) > 0:
            mod_redun = True
        if self.opt and not self.ts and not mod_redun and not self.no_freeze:
            input_line += "opt "
        elif self.opt and self.ts and not mod_redun and not self.no_freeze:
            input_line += "opt=(ts,calcfc,noeigentest) "
        elif self.opt and not self.ts and not mod_redun and self.no_freeze:
            input_line += "opt=NoFreeze "
        elif self.opt and self.ts and not mod_redun and self.no_freeze:
            input_line += "opt=(ts,calcfc,noeigentest,NoFreeze) "
        elif self.opt and not self.ts and mod_redun:
            input_line += "opt=ModRedundant "
        elif self.opt and self.ts and mod_redun:
            input_line += "opt=(ts,calcfc,noeigentest,ModRedundant) "
        if self.irc_opt:
            input_line += f"opt=(calcfc,maxstep=10,verytight) "
        if self.freq:
            input_line += "freq "
        input_line += f"{self.functional}/{self.basis} "
        if self.irc:
            input_line += f"irc=(lqa,rcfc,{self.irc_direction},maxcycle=50,maxpoint=100) guess=read "
        if self.dispersion_method:
            input_line += f"EmpiricalDispersion={self.dispersion_method} "
        if self.solvent:
            input_line += f"scrf=(SMD,solvent={self.solvent}) "
        if self.geom_chk:
            input_line += f"geom=AllCheck "
        if self.nbo:
            input_line += "pop=nboread "
        for option in self.route:
            input_line += f"{option} "
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
        if not self.geom_chk:
            xyz_directory = os.path.join(self.directory,self.xyz)
            with open(xyz_directory, "r") as f:
                xyz_lines = f.readlines()
                self.xyz_lines = xyz_lines[2:]
                f.close()
        else:
            pass