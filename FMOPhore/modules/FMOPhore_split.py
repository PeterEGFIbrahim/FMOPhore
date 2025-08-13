from Bio.PDB import *
import numpy as np
from Bio.PDB import PDBParser, Selection, Structure, Model, Chain, PDBIO, is_aa, Select
import sys
import os
import re
import pymol 
from pymol import cmd
import glob
import shutil
import Bio.PDB
import math
import argparse
from .FMOPhore_utility import EnvironmentGuard
from rdkit import Chem
from rdkit.Chem import AllChem
# EnvironmentGuard().enforce()

###################################################################################################
def distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2) 
###################################################################################################
#  SplitProcessor
###################################################################################################
class SplitProcessor:
    def __init__(self, pdb_file, distance_cutoff, same_target=False, Binding_Energy=False):
        self.pdb_file = pdb_file
        self.distance_cutoff = distance_cutoff if distance_cutoff != "no_cutoff" else np.inf
        self.same_target = same_target
        self.Binding_Energy = Binding_Energy
        self.chain_id_list = []
        self.ligands_list = []
    def process_FMO(self, ligand_info):
        if self.same_target:
            self.generate_cutoff_pdb_complex(distance_cutoff=self.distance_cutoff, ligand_info=ligand_info, same_target=True)
            self.process_FMO_complex(ligand_info)
        else:
            self.generate_cutoff_pdb_complex(distance_cutoff=self.distance_cutoff, ligand_info=ligand_info)
            self.process_FMO_complex(ligand_info)
        self.pymol_protein_add_H(self.distance_cutoff, ligand_info)
        self.modify_complex(self.distance_cutoff, ligand_info)
        # if self.Binding_Energy:
        self.process_FMO_protein(ligand_info)
        # os.chdir("..")
    def process_FMO_complex(self,ligand_info):
        pdb_lines = self.load_FMO_pdb_file(ligand_info)
        self.save_FMO_pdb_file(pdb_lines, ligand_info)
        self.renumber_FMO_complex_atoms(ligand_info)
    def process_FMO_protein(self,ligand_info):
        protein_lines = self.load_protein_file()
        self.save_outputprotein_file(protein_lines) 
        protein_lines = self.load_outputprotein_file()
        self.save_FMO_protein(protein_lines, ligand_info)
        self.renumber_FMO_protein_atoms(ligand_info)
        self.save_FMO_protein(protein_lines, ligand_info)
    def process_FMO_ligand(self,ligand_info):
        ligand_lines = self.load_ligand_file(ligand_info)
        self.save_FMO_ligand(ligand_lines, ligand_info)
    def split_complex(self):
        pdb_lines = self.load_pdb_file()
        ligands_list = [self.pdb_file[-13:-4]]
        for ligand in ligands_list:
            self.process_FMO(ligand)
    ################################################################
    def load_pdb_file(self):
        with open(self.pdb_file, 'r') as f:
            return f.readlines()
    ################################################################
    def generate_cutoff_pdb_complex(self, distance_cutoff, ligand_info, same_target=None):
        def _read_cofactor_name_from_log(log_path="../FMOPhore.log"):
            try:
                with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
                    for s in f:
                        s = s.strip()
                        if s.lower().startswith("cofactor:"):
                            # handles both "cofactor: SAH A 701" and "cofactor: SAH-A-701"
                            m = re.search(r'cofactor:\s*([A-Za-z0-9]{1,3})\b', s, re.I)
                            if m:
                                return m.group(1).upper()
            except FileNotFoundError:
                pass
            return None

        residues_names = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','HIE',
                          'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        _res_set = set(residues_names)

        _raw_cof = _read_cofactor_name_from_log("../FMOPhore.log")
        cofactor_name = _raw_cof if (_raw_cof and _raw_cof not in _res_set) else None

        if same_target:
            residue_file = '../../library_analysis/residues/residues.txt'
            residues_to_select = []
            with open(residue_file, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line:
                        residues_to_select.append(line)
        else:
            residue_file = f'{self.pdb_file[:-4]}_residues.txt'
            residues_to_select = []
            with open(residue_file, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line:
                        residues_to_select.append(line)
        pdb_file_path = os.path.join(self.pdb_file[:-4] + f'.pdb')
        parser = PDBParser(QUIET=True)
        renumbered_structure = parser.get_structure('output', pdb_file_path)
        ligand = None 
        ligand_chain = None
        ligand_numeber = None
        for model in renumbered_structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == "LIG":
                        ligand = residue
                        ligand_chain = chain.id
                        ligand_number = residue.id[1]
        residues = []
        selected_residues = set() 
        for model in renumbered_structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == ' ' or residue.get_id()[0] == 'H' or residue.get_resname() == 'HOH' or residue.get_resname() == cofactor_name:
                        for atom1 in ligand:
                            for atom2 in residue:
                                dist = atom1 - atom2 
                                if dist <= distance_cutoff:
                                    residue_number = residue.get_id()[1]
                                    residue_name = residue.get_resname()
                                    residue_chain_id = chain.id

                                    residue_key = (residue_name, residue_chain_id, residue_number)
                                    if residue_key not in selected_residues:
                                        residues.append((residue_name, residue_chain_id, residue_number, residue))
                                        selected_residues.add(residue_key)

                                        prev_resnum = residue_number - 1
                                        next_resnum = residue_number + 1
                                        for neighbor_resnum in [prev_resnum, next_resnum]:
                                            neighbor_key = (residue_name, residue_chain_id, neighbor_resnum)
                                            if neighbor_key not in selected_residues:
                                                try:
                                                    neighbor_residue = chain[neighbor_resnum]
                                                    residues.append((residue_name, residue_chain_id, neighbor_resnum, neighbor_residue))
                                                    selected_residues.add(neighbor_key)
                                                except KeyError:
                                                    pass
        if same_target:
            for residue_info in residues_to_select:
                residue_name, residue_chain_id, residue_number = residue_info.split('-')
                residue_number = int(residue_number)

                for model in renumbered_structure:
                    for chain in model:
                        if chain.get_id() == residue_chain_id:
                            try:
                                residue = chain[residue_number]
                                if residue.get_id()[0] == ' ' or residue.get_id()[0] == 'H' or residue.get_resname() == 'HOH' or residue.get_resname() == cofactor_name:
                                    found = True
                                    if residue_number not in selected_residues:
                                        residues.append((residue_name, residue_chain_id, residue_number, residue))
                                        selected_residues.add(residue_key)
                                    break
                            except KeyError:
                                pass
        else:
            pass
        def custom_sort_key(residue):
            chain_id, residue_number = residue[1], residue[2]
            return (chain_id, residue_number)
        residues.sort(key=custom_sort_key)
        sorted_residues = sorted(selected_residues, key=custom_sort_key)
        with open(f'{self.pdb_file[:-4]}_residues.txt', 'w') as file:
            for residue in sorted_residues:
                residue_name, chain_id, residue_number = residue 
                file.write(f"{residue_name}-{chain_id}-{residue_number}\n")
        with open(f'{self.pdb_file[:-4]}_residues.txt', 'r') as file:
            lines = file.read().splitlines()
            unique_lines = set()
            filtered_lines = []
            for line in lines:
                parts = line.split('-')
                chain_residue = f"{parts[1]}-{parts[2]}"
                if chain_residue not in unique_lines:
                    unique_lines.add(chain_residue)
                    filtered_lines.append(line)
        correct_residue_names = {}
        with open(self.pdb_file[:-4] + f'.pdb', 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith('ATOM'):
                    chain_id = line[21]
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    residue_key = f"{chain_id}-{res_num}"
                    correct_residue_names[residue_key] = res_name
        updated_lines = []
        for line in filtered_lines:
            parts = line.split('-')
            chain_residue = f"{parts[1]}-{parts[2]}"
            if chain_residue in correct_residue_names:
                parts[0] = correct_residue_names[chain_residue]
                updated_line = '-'.join(parts)
                updated_lines.append(updated_line)
        with open(f'{self.pdb_file[:-4]}_residues.txt', 'w') as file:
            file.write('\n'.join(updated_lines))
        unique_substrings = set()
        output_pdb_file = os.path.join(self.pdb_file[:-4] + f"_FMO_" + str(distance_cutoff) + ".pdb")
        with open(output_pdb_file, 'w') as pdb_file:
            for residue_tuple in residues:
                residue_object = residue_tuple[3] 
                chain_id = residue_tuple[1]  
                residue_number = residue_tuple[2]  
                for atom in residue_object:
                    atom_name = atom.name
                    atom_serial_number = atom.serial_number
                    atom_coordinates = atom.coord
                    pdb_line = f"ATOM  {atom_serial_number:5d} {atom_name:<4} {residue_object.resname:3} {chain_id:1}{residue_number:4d}    {atom_coordinates[0]:8.3f}{atom_coordinates[1]:8.3f}{atom_coordinates[2]:8.3f}{1.00:6.2f}{0.00:6.2f}          {atom_name[0]:>2}\n"
                    pdb_file.write(pdb_line)
            for atom in ligand:
                atom_name = atom.name
                atom_serial_number = atom.serial_number
                atom_coordinates = atom.coord
                pdb_line = f"HETATM{atom_serial_number:5d} {atom_name:<4} {ligand.resname:3} {ligand_chain:1}{ligand_number:4d}    {atom_coordinates[0]:8.3f}{atom_coordinates[1]:8.3f}{atom_coordinates[2]:8.3f}{1.00:6.2f}{0.00:6.2f}          {atom_name[0]:>2}\n"
                pdb_file.write(pdb_line)
        def correct_HETATM(input_file, output_file):
            hetatm_lines = []
            other_lines = []
            with open(input_file, 'r') as infile:
                for line in infile:
                    resname = line[17:20].strip().upper()
                    if resname == "HOH" or resname == cofactor_name:
                        line = line[:0] + "HETATM" + line[6:]
                    if line.startswith("HETATM"):
                        hetatm_lines.append(line)
                    else:
                        if "END" not in line:
                            other_lines.append(line)
            with open(output_file, 'w') as outfile:
                outfile.writelines(other_lines)
                outfile.writelines(hetatm_lines)
            renumbered_lines = []
            with open(input_file, 'r') as infile:
                line_counter = 0
                for line in infile:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        line_counter += 1
                        line = line[:6] + str(line_counter).rjust(5) + line[11:]
                        renumbered_lines.append(line)
            with open(output_file, 'w') as outfile:
                outfile.writelines(renumbered_lines)
        correct_HETATM(output_pdb_file, output_pdb_file)
    ################################################################
    def load_FMO_pdb_file(self, ligand_info):
        output_pdb_file = self.pdb_file[:-4] + f"_FMO_" + str(self.distance_cutoff) + ".pdb"
        with open(output_pdb_file, 'r') as f:
            return f.readlines()
    ################################################################
    def save_FMO_pdb_file(self, pdb_lines, ligand_info):
        output_pdb_file = self.pdb_file[:-4] + f"_FMO_" + str(self.distance_cutoff) + ".pdb"
        with open(output_pdb_file, 'w') as f:
            f.writelines(pdb_lines)
    ################################################################
    def renumber_FMO_complex_atoms(self, ligand_info):
        parser = PDBParser(QUIET=True)
        pdb_file_path = os.path.join(self.pdb_file[:-4] + f"_FMO_" + str(self.distance_cutoff) + ".pdb")
        renumbered_structure = parser.get_structure('renumbered', self.pdb_file[:-4] + f"_FMO_" + str(self.distance_cutoff) + ".pdb")
        output_filename = self.pdb_file[:-4] + f"_FMO_" + str(self.distance_cutoff) + ".pdb"
        with open(pdb_file_path, 'w') as f:
            writer = PDBIO()
            writer.set_structure(renumbered_structure)
            writer.save(f)
    ################################################################
    def renumber_FMO_protein_atoms(self, ligand_info):
        parser = PDBParser(QUIET=True)
        pdb_file_path = os.path.join(self.pdb_file[:-4] + f"_FMO_" + "protein.pdb")
        renumbered_structure = parser.get_structure('renumbered', self.pdb_file[:-4] + f"_FMO_" + "protein.pdb")
        output_filename = self.pdb_file[:-4] + f"_FMO_" + "protein.pdb"
        with open(pdb_file_path, 'w') as f:
            writer = PDBIO()
            writer.set_structure(renumbered_structure)
            writer.save(f)
    ################################################################
    # def pymol_protein_add_H(self, distance_cutoff, ligand_info): 
    #     folder_path = '.' 
    #     original_PDB_file = glob.glob(self.pdb_file[:-4] + f"_FMO_" + str(self.distance_cutoff) + ".pdb")
    #     for file in original_PDB_file:
    #         new_filename = file[:-4] + '_H.pdb'
    #         shutil.copy(file, new_filename)
    #         ligand_line = next(line for line in open(new_filename) if line[17:20].strip() == 'LIG')
    #         cmd.load(new_filename)
    #         cmd.set_name(new_filename.split('.')[0], "bla")
    #         cmd.select('mysele', 'name C')
    #         cmd.h_add("mysele")
    #         # cmd.select('mysele', 'name CB')
    #         # cmd.h_add("mysele")
    #         cmd.select('mysele', 'name N')
    #         cmd.h_add('mysele')
    #         cmd.select('mysele', 'resn LIG and not name CL* and not name N* and not name O*')
    #         cmd.h_add('mysele')
    #         cmd.save(new_filename)

    def pymol_protein_add_H(self, distance_cutoff, ligand_info): 
        original_pdb_files = glob.glob(self.pdb_file[:-4] + f"_FMO_{distance_cutoff}.pdb")
        for file in original_pdb_files:
            new_filename = file[:-4] + '_H.pdb'
            shutil.copy(file, new_filename)
            # Load structure into PyMOL
            cmd.reinitialize() # Clear any previous session
            cmd.load(new_filename, "mol")
            # Only select ATOM records (i.e., standard protein atoms)
            # Select backbone C and N atoms from protein (ATOM entries only)
            cmd.select("protein_c", "model mol and polymer.protein and name C")
            cmd.h_add("protein_c")
            cmd.select("protein_n", "model mol and polymer.protein and name N")
            cmd.h_add("protein_n")
            cmd.save(new_filename, "mol")

    # def pymol_protein_add_H(self, distance_cutoff, ligand_info):
    #     # Locate the original file
    #     input_files = glob.glob(self.pdb_file[:-4] + f"_FMO_{distance_cutoff}.pdb")
    #     if not input_files:
    #         raise FileNotFoundError("Input PDB file not found.")
    #     for input_file in input_files:
    #         output_file = input_file[:-4] + '_H.pdb'
    #         # Copy the original file for safety (optional)
    #         shutil.copy(input_file, output_file)
    #         # Load without sanitization to preserve coordinates
    #         mol = Chem.MolFromPDBFile(output_file, removeHs=False, sanitize=False)
    #         if mol is None:
    #             raise ValueError(f"Could not parse the PDB file: {output_file}")
    #         # Add hydrogens, preserving original coordinates
    #         mol_with_H = Chem.AddHs(mol, addCoords=True)
    #         # ✅ Do NOT call AllChem.EmbedMolecule — it will override coordinates
    #         # ✅ Coordinates are already present, just extend to added Hs
    #         # Save the result
    #         Chem.MolToPDBFile(mol_with_H, output_file)
    #         print(f"Hydrogens added with original coordinates preserved: {output_file}")

    ################################################################
    def modify_complex(self, distance_cutoff, ligand_info):
        f2_mod_filename = self.pdb_file[:-4] + f"_FMO_" + str(self.distance_cutoff) + "_H.pdb"    
        with open(f2_mod_filename, "r") as file:
            lines = file.readlines()
        amino_acid_mapping = {
            " NZ  LYS": "N1+",
            " O   ASP": "O1-",
            " N   ARG": "N1+",
            # " NH1 ARG": "N1+",
            # " NE  ARG": "N1+",
            " O   GLU": "O1-"
        }
        his_line_count = {}
        atom_names_by_residue = {}
        for line in lines:
            if "HIS" in line[17:20]:
                chain_id = line[20:21]
                residue_id = line[21:26]
                atom_name = line[12:16].strip()
                if residue_id not in atom_names_by_residue:
                    atom_names_by_residue[residue_id] = []
                atom_names_by_residue[residue_id].append(atom_name)  
        his_line_count = {}
        required_atoms = {'HD1', 'HD2','HE1', 'HE2', 'HD4', 'HD5','HE4', 'HE5'}
        for residue_id, atoms in atom_names_by_residue.items():
            present_atoms_count = len(set(atoms) & required_atoms)
            if present_atoms_count >= 4:
                his_line_count[residue_id] = his_line_count.get(residue_id, 0) + 1
                amino_acid_mapping[f" N   HIS{chain_id}{residue_id}"] = "N1+"
        final_modified_lines = []
        for line in lines:
            if line[12:20] in amino_acid_mapping:
                modified_line = line[:77] + amino_acid_mapping[line[12:20]] + line[80:]
                final_modified_lines.append(modified_line)
            elif line[12:26] in amino_acid_mapping:
                modified_line = line[:77] + amino_acid_mapping[line[12:26]] + line[80:]
                final_modified_lines.append(modified_line)
            else:
                final_modified_lines.append(line)
        with open(f2_mod_filename, "w") as file:
            file.writelines(final_modified_lines)
        def rearrange_atom_file(input_file, file_path):
            atom_lines = []
            hoh_lines = []
            hetatm_lines = []
            other_lines = []
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('ATOM'):
                        if 'HOH' in line.split()[3]: 
                            hoh_lines.append(line)
                        else:
                            atom_lines.append(line)
                    elif line.startswith('HETATM'):
                        hetatm_lines.append(line)
                    else:
                        other_lines.append(line)
            with open(file_path, 'w') as file:
                file.writelines(atom_lines + hoh_lines + hetatm_lines + other_lines)
        def correct_ATOM(input_file, output_file):
            hetatm_lines = []
            other_lines = []
            with open(input_file, 'r') as infile:
                for line in infile:
                    if line[17:20] != "LIG" and line[0:6] == "HETATM":
                        line = "ATOM  " + line[6:]
                    other_lines.append(line)
            with open(output_file, 'w') as outfile:
                outfile.writelines(other_lines)
        correct_ATOM(f2_mod_filename,f2_mod_filename)
        def correct_HETATM(input_file, output_file):
            hetatm_lines = []
            other_lines = []
            with open(input_file, 'r') as infile:
                for line in infile:
                    if line[17:20] == "LIG":
                        line = line[:0] + "HETATM" + line[6:]
                    if line.startswith("HETATM"):
                        hetatm_lines.append(line)
                    else:
                        if "END" not in line:
                            other_lines.append(line+ "\n")
            with open(output_file, 'w') as outfile:
                outfile.writelines(other_lines)
                outfile.writelines(hetatm_lines)
            renumbered_lines = []
            with open(input_file, 'r') as infile:
                line_counter = 0
                for line in infile:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        line_counter += 1
                        line = line[:6] + str(line_counter).rjust(5) + line[11:]
                        renumbered_lines.append(line)
            with open(output_file, 'w') as outfile:
                outfile.writelines(renumbered_lines)
        correct_HETATM(f2_mod_filename,f2_mod_filename)
    ################################################################
    #  Extract the Protein
    ################################################################
        cwd = os.getcwd()
        # if self.Binding_Energy:
        class WaterAndProteinSelect(Select):
            def accept_residue(self, residue):
                ions_and_waters = ["MG", "ZN", "HOH"]
                return is_aa(residue) or residue.get_resname() in ions_and_waters
        def extract_protein(path, ligand_info):
            protein = None
            for pdb_file in os.listdir(path + "/"):
                if pdb_file.endswith(self.pdb_file[:-4] + f"_FMO_" + str(self.distance_cutoff) + "_H.pdb") and not pdb_file.startswith("lig_"):
                    pdb_code = pdb_file[:-4]
                    pdb = PDBParser().get_structure(pdb_code, path + "/" + pdb_file)
                    if protein is None:
                        protein = pdb.copy()
                    else:
                        for model in pdb:
                            protein.add(model)
            if protein is not None:
                io = PDBIO()
                io.set_structure(protein)
                io.save("protein.pdb", WaterAndProteinSelect())
    ################################################################
    #  Extract the Ligand
    ################################################################
        def is_het(residue):
            res = residue.id[0]
            return res != " " and res != "W"
        class ResidueSelect(Select):
            def __init__(self, chain, residue):
                self.chain = chain
                self.residue = residue
            def accept_chain(self, chain):
                return chain.id == self.chain.id
            def accept_residue(self, residue):
                return residue == self.residue and is_het(residue) 
        def extract_ligands(path, ligand_info):
            for pfb_file in os.listdir(path + "/"):
                if pfb_file.endswith(self.pdb_file[:-4] + f"_FMO_" + str(self.distance_cutoff) + "_H.pdb" ) and not pfb_file.startswith("lig_"):
                    pdb_code = pfb_file[:-4]
                    pdb = PDBParser().get_structure(pdb_code, path + "/" + pfb_file)
                    io = PDBIO()
                    io.set_structure(pdb)
                    for model in pdb:
                        for chain in model:
                            chain_id = chain.id
                            chain_ligands = [residue for residue in chain if is_het(residue)]
                            if not chain_ligands:
                                continue
                            for ligand in chain_ligands:
                                ligand_id = ligand.get_id()[0]
                                filename = f"{self.pdb_file[:-4]}_FMO_lig_H.pdb"
                                io.save(filename, ResidueSelect(chain, ligand))
        path = cwd
        extract_protein(path, ligand_info) 
        extract_ligands(path, ligand_info)
    ################################################################
    #  process the Protein
    ################################################################
    def load_protein_file(self):
        with open("protein.pdb", 'r') as f:
            return f.readlines()
    def save_outputprotein_file(self, protein_lines):
        output_file = 'output.pdb'  
        with open(output_file, 'w') as f:
            f.writelines(protein_lines)
    def load_outputprotein_file(self):
        os.remove('protein.pdb')
        with open("output.pdb", 'r') as f:
            return f.readlines() 
    def save_FMO_protein(self, protein_lines, ligand_info):
        output_file = 'atom_modified.pdb'
        if os.path.exists(output_file):
            os.remove(output_file)
        with open(output_file, 'w') as f:
            f.writelines(protein_lines)
        with open(output_file, "r") as file:
            lines = file.readlines()
        amino_acid_mapping = {
            " NZ  LYS": "N1+",
            " O   ASP": "O1-",
            " N   ARG": "N1+",
            # " NH1 ARG": "N1+",
            # " NE  ARG": "N1+",
            " O   GLU": "O1-"
        }
        his_line_count = {}
        atom_names_by_residue = {}
        for line in lines:
            if "HIS" in line[17:20]:
                chain_id = line[20:21]
                residue_id = line[21:26]
                atom_name = line[12:16].strip()
                if residue_id not in atom_names_by_residue:
                    atom_names_by_residue[residue_id] = []
                atom_names_by_residue[residue_id].append(atom_name)  
        his_line_count = {}
        required_atoms = {'HD1', 'HD2','HE1', 'HE2', 'HD4', 'HD5','HE4', 'HE5'}
        for residue_id, atoms in atom_names_by_residue.items():
            present_atoms_count = len(set(atoms) & required_atoms)
            if present_atoms_count >= 4:
                his_line_count[residue_id] = his_line_count.get(residue_id, 0) + 1
                amino_acid_mapping[f" N   HIS{chain_id}{residue_id}"] = "N1+"
        final_modified_lines = []
        for line in lines:
            if line[12:20] in amino_acid_mapping:
                modified_line = line[:77] + amino_acid_mapping[line[12:20]] + line[80:]
                final_modified_lines.append(modified_line)
            elif line[12:26] in amino_acid_mapping:
                modified_line = line[:77] + amino_acid_mapping[line[12:26]] + line[80:]
                final_modified_lines.append(modified_line)
            else:
                final_modified_lines.append(line)
        output_filename = f"{self.pdb_file[:-4]}_FMO_protein.pdb"
        with open(output_filename, "w") as file:
            file.writelines(final_modified_lines)
    ################################################################
    #  process the Ligand
    ################################################################
    def load_ligand_file(self, ligand_info):
        ligand_file = next((filename for filename in os.listdir() if "lig" in filename), None)
        if ligand_file is not None:
            with open(ligand_file, 'r') as f:
                return f.readlines()
        else:
            print("No ligand file found.")
            return []
    def save_FMO_ligand(self, ligand_lines, ligand_info):
        with open(f"../" + f"{self.pdb_file[:-14]}.pdb", "r") as f1:
            f1_lines = f1.readlines()
        new_f2_lines = []
        for f2_line in ligand_lines:
            f2_prefix = f2_line[32:67]
            f1_matching_line = next((f1_line for f1_line in f1_lines if f1_line[32:67] == f2_prefix), None)
            if f1_matching_line is not None:
                f1_suffix = f1_matching_line[67:]
                f2_suffix = f2_line[67:]
                new_f2_line = f2_line[:32] + f2_prefix + f1_suffix
            else:
                new_f2_line = f2_line
            new_f2_lines.append(new_f2_line)
        with open(f"{self.pdb_file[:-4]}_FMO_lig_H.pdb", 'w') as f:
            f.writelines(new_f2_lines)
###################################################################################################
if __name__ == "__main__":
    """
    FMOPhore V.0.1 - SplitProcessor - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    LIG_processor = SplitProcessor(args.pdb_file)
    LIG_processor.split_complex()
