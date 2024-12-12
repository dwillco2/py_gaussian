import os
from .config import *

class GenJob:
    def __init__(self, 
                 directory, 
                 name,
                 xyz,
                 charge,
                 multiplicity,
                 nprocs, 
                 memory,
                 solvent,
                 ) -> None:
        self.directory = directory
        self.name = name
        if not xyz:
            self.xyz = f"{self.name}.xyz"
        else:
            self.xyz = xyz
        self.charge = charge
        self.multiplicity = multiplicity
        self.nprocs = nprocs
        self.memory = memory
        self.solvent = solvent
        if not self.charge:
            self.get_charge()
        self.read_xyz()

    def gen_output(self):
        pass

    def gen_bash(self, runtime=48):
        pass

    def gen_bash_template(self, runtime):
        template = bash_template.replace("_JOBNAME_",self.output_name)
        template = template.replace("_RUNTIME_",f"{runtime}")
        template = template.replace("_MEMORY_",f"{self.memory}")
        template = template.replace("_NPROCS_",f"{self.nprocs}")
        return template
    
    
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
            
    
    def read_xyz(self):
        if self.xyz:
            xyz_directory = os.path.join(self.directory,self.xyz)
            with open(xyz_directory, "r") as f:
                xyz_lines = f.readlines()
                self.xyz_lines = xyz_lines[2:]
                f.close()