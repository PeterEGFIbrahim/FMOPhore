import os
import glob
import shutil
import argparse
from Bio import PDB
from Bio.PDB import Superimposer
from .FMOPhore_utility import EnvironmentGuard
# EnvironmentGuard().enforce()

def rename_pdb_files(directory_path):
    pdb_files = [f for f in os.listdir(directory_path) if f.endswith(".pdb")]
    for pdb_file in pdb_files:
        file_path = os.path.join(directory_path, pdb_file)
        with open(file_path, 'rb') as f:
            lines = f.readlines()
        updated_lines = []
        for line in lines:
            line = line.decode('utf-8', errors='ignore')
            if line.startswith("HETATM"):
                updated_line = line.replace("UNK", "UNK")
                updated_lines.append(updated_line.encode('utf-8'))
            else:
                updated_lines.append(line.encode('utf-8'))
        with open(file_path, 'wb') as f:
            f.writelines(updated_lines)
def align_pdb_backbone(reference_pdb_file, pdb_file, output_file):
    parser = PDB.PDBParser(QUIET=True)
    structure1 = parser.get_structure("reference", reference_pdb_file)
    structure2 = parser.get_structure("structure", pdb_file)
    atoms1 = []
    atoms2 = []
    for model1, model2 in zip(structure1, structure2):
        for chain1, chain2 in zip(model1, model2):
            for residue1, residue2 in zip(chain1, chain2):
                if PDB.is_aa(residue1) and PDB.is_aa(residue2):
                    for atom_name in ["N", "CA", "C"]:
                        atom1 = residue1[atom_name]
                        atom2 = residue2[atom_name]
                        atoms1.append(atom1)
                        atoms2.append(atom2)
    superimposer = Superimposer()
    superimposer.set_atoms(atoms1, atoms2)
    superimposer.apply(structure2)
    io = PDB.PDBIO()
    io.set_structure(structure2)
    io.save(output_file)
def add_charges(before_aligned, aligned):
    with open(before_aligned, "r") as f1:
        f1_lines = f1.readlines()
    new_f2_lines = []
    with open(aligned, "r") as f2:
        f2_lines = f2.readlines()
    for f2_line in f2_lines:
        f2_prefix = f2_line[12:26]
        f1_matching_line = next((line for line in f1_lines if line[12:26] == f2_prefix), None)
        if f1_matching_line:
            f1_suffix = f1_matching_line[67:]
            f2_coords = f2_line[26:67]
            new_f2_line = f2_line[:67] + f1_suffix
        else:
            new_f2_line = f2_line
        new_f2_lines.append(new_f2_line)
    with open(aligned, "w") as f2_mod:
        f2_mod.writelines(new_f2_lines)
def remove_spaces_between_lines_and_END(input_file, output_file):
    hetatm_lines = []
    other_lines = []
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith("HETATM"):
                hetatm_lines.append(line)
            else:
                if "END" not in line:
                    other_lines.append(line)
    with open(output_file, 'w') as outfile:
        outfile.writelines(other_lines)
        outfile.writelines(hetatm_lines)
    with open(output_file, 'r') as f:
        lines = f.readlines()
    cleaned_lines = [line.strip() for line in lines if line.strip()]
    with open(output_file, 'w') as f:
        f.write('\n'.join(cleaned_lines))
def aligner(reference_pdb_file, pdb_files_pattern, output_dir):
    """
    FMOPhore V.0.1 - Aligner - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if os.path.isdir(pdb_files_pattern):
        pdb_files = glob.glob(os.path.join(pdb_files_pattern, "*.pdb"))
    else:
        pdb_files = glob.glob(pdb_files_pattern)
    for input_path in pdb_files:
        base_filename = os.path.basename(input_path)
        corresponding_output_path = os.path.join(output_dir, base_filename)
        align_pdb_backbone(reference_pdb_file, input_path, corresponding_output_path)
        add_charges(input_path, corresponding_output_path)
        remove_spaces_between_lines_and_END(corresponding_output_path, corresponding_output_path)
        os.remove(input_path)
        shutil.move(corresponding_output_path, input_path)
    shutil.rmtree(output_dir)
