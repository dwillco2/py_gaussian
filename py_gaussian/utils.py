from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
import os
import numpy as np
import pandas as pd

def embed_and_optimize(smiles):
    # generate UFF optimised mol from smiles
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    if mol:
        rdDistGeom.EmbedMolecule(mol, params=AllChem.ETKDGv3())
        AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
        return mol


def save_xyz_and_charge(mol, path):
    # save xyz and charge from optimised mol
    _, tail = os.path.split(path)
    new_path = os.path.join(path,tail)
    os.mkdir(path)
    Chem.MolToXYZFile(mol,f"{new_path}.xyz")
    charge = Chem.GetFormalCharge(mol)
    with open(f"{new_path}.CHRG", "w") as f:
        f.write(f"{charge}")
        f.close


def is_optimised(file):
    # returns True if gaussian .log file geomertry optimisation is complete
    with open(file, "r") as f:
        lines = f.readlines()
        f.close()
    optimisation_complete = False
    stationary_point = False
    imaginary_frequencies = False
    ts = False
    n_imaginary_frequencies = 0
    frequencies_complete = False
    for line in lines:
        if "Search for a saddle point of order  1." in line:
            ts = True
        elif "Optimization completed" in line:
            optimisation_complete = True
        elif "Stationary point found." in line:
            stationary_point = True
        elif "Frequencies" in line:
            split_line = line.split()
            for i in split_line[2:]:
                if float(i) < 0:
                    imaginary_frequencies = True
                    n_imaginary_frequencies += 1
        elif "Harmonic frequencies" in line:
            frequencies_complete = True
    
    if optimisation_complete and stationary_point and not imaginary_frequencies and frequencies_complete and not ts:
        return True
    elif optimisation_complete and stationary_point and imaginary_frequencies and ts and n_imaginary_frequencies == 1 and frequencies_complete:
        return True
    else:
        return False


def freq_finished(file):
    # returns True if Gaussian .log file frequency optimisation is complete
    with open(file, "r") as f:
        lines = f.readlines()
        f.close()
    frequencies_complete = False
    for line in lines:
        if "Harmonic frequencies" in line:
            frequencies_complete = True
    if frequencies_complete:
        return True
    else:
        return False
    

def read_gaussian_output(file):
    # takes a Gaussian .log file and reads SCF and thermochemistry info into dataframe
    file_name = os.path.split(file)
    name = file_name[1].replace(".log","")
    with open(file, "r") as f:
        lines = f.readlines()
        f.close()
    scf = 0
    zpe = 0
    e = 0
    h = 0
    g = 0
    failed = False

    for line in lines:
        if "SCF Done:" in line:
            split_scf = line.split()
            scf = float(split_scf[4])
        elif "Sum of electronic and zero-point Energies=" in line:
            zpe_split = line.split()
            zpe = float(zpe_split[6])
        elif "Sum of electronic and thermal Energies=" in line:
            e_split = line.split()
            e = float(e_split[6])
        elif "Sum of electronic and thermal Enthalpies=" in line:
            h_split = line.split()
            h = float(h_split[6])
        elif "Sum of electronic and thermal Free Energies=" in line:
            g_split = line.split()
            g = float(g_split[7])
        elif "Error termination via" in line:
            failed = True
        else:
            pass
    
    df = pd.DataFrame(data={"name": [name],"scf": [scf], "zpe": [zpe], "E": [e], "H": [h], "G": [g], "Failed": [failed]})
    return df
    

def split_xyz(file):
    # takes multi-xyz as input and splits them into seperate xyz lists
    with open(file,"r") as f:
        lines = f.readlines()
        f.close()
    split_lines = []
    xyz_n = -1
    first_line = False
    second_line = False
    for i,line in enumerate(lines):
        if line.split()[0].isdigit():
            split_lines.append([])
            xyz_n += 1
            split_lines[xyz_n].append(line)
            first_line = True
        elif first_line:
            split_lines[xyz_n].append(line)
            first_line = False
            second_line = True
        elif second_line:
            split_lines[xyz_n].append(line)
    return split_lines


