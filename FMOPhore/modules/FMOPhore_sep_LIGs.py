import argparse
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
import subprocess
import csv
import math
from pathlib import Path
import multiprocessing
###################################################################################################
#  LIGProcessor
###################################################################################################
class LIGProcessor:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
    def chain_id_list_all(self):
        return self.chain_id_list
    def ligands_list_all(self):
        return self.ligands_list
    ################################################################
    chain_id_list = []
    ligands_list  = []
    def process_ligands(self, ligand):
        ligand_info = "_".join(ligand.split())
        pdb_lines_ligand = self.load_processed_pdb_file(ligand_info)
        pdb_lines_ligand = self.filter_lines(pdb_lines_ligand)
        pdb_lines_ligand = self.include_ions(pdb_lines_ligand)
        pdb_lines_ligand = self.move_conect_lines_to_end(pdb_lines_ligand)
        pdb_lines_ligand = self.move_END_lines_to_end(pdb_lines_ligand)
        self.save_processed_pdb_file(pdb_lines_ligand, ligand_info)
        def check_multiple_chains_same_ligand(log_path):
            with open(log_path, 'r') as file:
                for line in file:
                    if " falls between multiple chains " in line:
                        return True
            return False
        log_path = "FMOPhore.log"
        if check_multiple_chains_same_ligand(log_path):
            self.sep_chains(ligand_info, multiple_chains=True)
        else:
            self.sep_chains(ligand_info, multiple_chains=False)
        ligand_info_folder = f"{self.pdb_file[:-4]}_{ligand_info}"
        with os.scandir(ligand_info_folder):
            os.chdir(ligand_info_folder)
            pdb_lines_ligand = self.load_processed_pdb_file(ligand_info)
            pdb_lines_ligand = self.rename_LIG(pdb_lines_ligand)
            pdb_lines_ligand = self.modify_atom_names(pdb_lines_ligand)
            self.save_processed_pdb_file_ligand_info(pdb_lines_ligand, ligand_info)
            self.renumber_important_residues(ligand_info)
            print(True)
        os.chdir("..")
    def process_complex_LIG(self):
        pdb_lines = self.load_pdb_file()
        pdb_lines = self.separate_ligands(pdb_lines)
        for ligand in set(self.ligands_list):
            self.process_ligands(ligand)
    ################################################################
    def load_pdb_file(self):
        with open(self.pdb_file, 'r') as f:
            return f.readlines()
    ################################################################
    def separate_ligands(self, pdb_lines):
        excluded_strings = [
          'HD2 ASP','FMT', 'ANISOU', 'BU3','NAG' , 'PO4',
          'PEG', 'GOL', 'BO3', 'EDO', 'SO4',  'NO3',
          'ACE', 'NMA', 'NME', 'ACT', 'MES', 'TRS',
          'OCY', "HA2 PHE", 'SEP', 'BME', 'CSO', 'IMD',
          'TPO', 'TCE', 'MDP', 'IOD', 'NI', 'IPA', 'CO',
          'ZN', 'CL','T3P', 'WAT', 'HOH', 'MG'
          ]
        for line in pdb_lines:
            record_type = line[:6].strip()
            if record_type == 'HETATM' and line[17:20].strip() not in excluded_strings:
                chain_id = line[21:22]
                ligand = line[17:26]
                self.ligands_list.append(ligand)
        ligands_by_chain = {}
        with open("FMOPhore.log", 'r') as f_in:
            for line in f_in:
                line = line.strip()
                if "falls in chain" in line:
                    ligand_name_from_log = line.split()[0]  
                    if ligand_name_from_log not in excluded_strings:
                        ligand_info = line[0:9]
                        chain_id = line[4:5]  
                        if chain_id not in ligands_by_chain:
                            ligands_by_chain[chain_id] = []
                        ligands_by_chain[chain_id].append(f"{ligand_info}")
        filtered_pdb_lines = pdb_lines[:]
        for ligands in ligands_by_chain.values():
            for ligand in ligands:
                filtered_pdb_lines = [line for line in filtered_pdb_lines if ligand not in line]
        if any(len(ligands) >= 1 for ligands in ligands_by_chain.values()):
            for chain_id, ligands in ligands_by_chain.items():
                for ligand in ligands:
                    ligand_info = "_".join(ligand.split())
                    pdb_file_name = f"{self.pdb_file[:-4]}_{ligand_info}.pdb"
                    with open(pdb_file_name, 'w') as f:
                        written_ligand = False
                        for line in filtered_pdb_lines:
                            f.write(line.rstrip() + '\n')
                        for line in pdb_lines:
                            if re.search(fr"\b{ligand}\b", line):
                                f.write(line.rstrip() + '\n')
                                written_ligand = True
    ################################################################
    def load_processed_pdb_file(self, ligand_info):
        pdb_file_path = f"{self.pdb_file[:-4]}_{ligand_info}.pdb"
        with open(pdb_file_path, 'r') as f:
            return f.readlines()
    ###############################################################
    def filter_lines(self, pdb_lines):
        excluded_strings = [
          'HD2 ASP','FMT', 'ANISOU',  'NO3','NAG' , 'PO4',
          'PEG', 'GOL', 'BO3', 'EDO', 'SO4', 'IPA', 
          'ACE', 'NMA', 'NME', 'ACT', 'MES', 'TRS', 'MLA',
          'OCY', "HA2 PHE", 'SEP', 'BME', 'CSO', 'IMD', 'CO', 
          'TPO', 'TCE', 'MDP', 'IOD', 'NI'
                     , 'SAH', 'MTA', 'MG'
          ]
        pdb_lines = [
            line for line in pdb_lines
            if not (any(substr in line for substr in excluded_strings) or
                line[17:20] in ['CD ',' CD', ' CL', 'CL ', 'NA ', ' NA'] or
                # line[17:20] in [' CD', 'NA '] or
                line[16:17] in ['B', 'C'] or
                line[26:27] == 'A' or 
                (line.startswith("ATOM") and line[12:16].strip() == 'H3' and line[17:20].strip() == 'THR') or
                (line[0:4] == "ATOM" and line[12:16].strip() == 'OXT'))]
        for i, line in enumerate(pdb_lines):
            record_type = line[:6].strip()
            residue_name = line[17:20].strip()
            if record_type == 'HETATM' and residue_name != 'HOH':
                if line[22].isdigit() or line[22].isalpha():
                    pdb_lines[i] = line[:22] + ' ' + line[23:]
        return pdb_lines
    ################################################################
    def include_ions(self, pdb_lines):
        def is_float(value):
            try:
                float(value)
                return True
            except ValueError:
                return False
        include_strings = ['MG', 'ZN', 'CO', 'CL']  
        for i, line in enumerate(pdb_lines):
            if any(substr in line[17:20].strip() for substr in include_strings) and is_float(line[30:38].strip()):
                pdb_lines[i] = "ATOM  " + line[6:]  
        return pdb_lines
    ################################################################
    def move_conect_lines_to_end(self, pdb_lines):
        conect_lines = [line for line in pdb_lines if line.startswith("CONECT")]
        pdb_lines = [line for line in pdb_lines if not line.startswith("CONECT")]
        pdb_lines.extend(conect_lines)
        return pdb_lines
    def move_END_lines_to_end(self, pdb_lines):
        conect_lines = [line for line in pdb_lines if line.startswith("END")]
        pdb_lines = [line for line in pdb_lines if not line.startswith("END")]
        pdb_lines.extend(conect_lines)
        return pdb_lines
    ################################################################
    def save_processed_pdb_file(self, pdb_lines, ligand_info):
        with open(f"{self.pdb_file[:-4]}_{ligand_info}.pdb", 'w') as f:
            f.writelines(pdb_lines)
    ################################################################
    def correct_ligand_chain(self, pdb_lines):
        excluded_strings = [
          'WAT', 'T3P',  'FMT', 'ANISOU',  'NO3','NAG' , 'PO4',
          'PEG', 'GOL', 'BO3', 'EDO', 'SO4', 'ACE', 'TRS',
          'NMA', 'NME', 'ACT', 'MES', 'OCY','SEP', 'HOH',
          'TPO', 'TCE', 'MDP', 'IOD', 'MG'
                     , 'SAH', 'MTA'
        ]
        ligands = {}
        for line in pdb_lines:
            if line.startswith("HETATM") and line[17:20].strip() not in excluded_strings:
                ligand_key = (line[17:20], line[22:26])
                chain_id = line[21]
                if ligand_key in ligands:
                    all_chains = [chain for chains in ligands.values() for chain in chains]
                    most_common_chain = max(set(all_chains), key=all_chains.count)
                else:
                    ligands[ligand_key] = chain_id
        with open("FMOPhore.log", "a") as error_file:  
            ligand_dict = {}
            for i, line in enumerate(pdb_lines):
                if line.startswith("HETATM"):
                    ligand = line[17:26]
                    ligand_name = line[17:20]
                    ligand_num = line[22:26]
                    chain_id = line[21]
                    ligand_key = (ligand_name, ligand_num)
                    if ligand_key in ligand_dict:
                        existing_chain_id = ligand_dict[ligand_key]
                        if chain_id != existing_chain_id:
                            pdb_lines[i] = line[:21] + most_common_chain + line[22:]
                            if ligand_key[0].upper() not in [' ZN', 'T3P', 'WAT', 'HOH']:
                                error_message = "{} falls between multiple chains {} !".format(ligand, chain_id)
                                error_file.write(error_message + "\n")
                        else:
                            if ligand_key[0].upper() not in [' ZN', 'T3P', 'WAT', 'HOH']:
                                error_message = "{} falls in chain {} !".format(ligand, chain_id)
                                error_file.write(error_message + "\n")
                    else:
                        ligand_dict[ligand_key] = chain_id
        seen_lines = set()
        filtered_lines = []
        with open("FMOPhore.log", 'r') as f_in:
            for line in f_in:
                line = line.strip()
                if line not in seen_lines:
                    filtered_lines.append(line)
                    seen_lines.add(line)
        with open("FMOPhore.log", 'w') as f_out:
            for line in filtered_lines:
                f_out.write(line + '\n')
        ligands_count = 0
        ligand_chain_combinations = set()
        for line in pdb_lines:
            record_type = line[:6].strip()
            if record_type == 'HETATM' and line[17:20].strip() not in [' ZN', 'T3P', 'WAT', 'HOH']:
                chain_id = line[21:22]
                ligand = line[17:26]
                self.chain_id_list.append(chain_id)  
                self.ligands_list.append(ligand)
        return pdb_lines
    ################################################################
    def ligand_charge(self, pdb_lines):
        with open("FMOPhore.log", 'r') as log_file:
            log_lines = log_file.readlines()
        charge_dict = {'O1+': '+1', 'N1+': '+1', 'O1-': '-1', 'S1-': '-1', 'l1-': '-1', 'N1-': '-1'}
        ligand_charges = {} 
        for i, line in enumerate(log_lines):
            ligand_name = line[:9] 
            if ligand_name not in ligand_charges:
                ligand_charges[ligand_name] = 0 
            for pdb_line in pdb_lines:
                if pdb_line.startswith('HETATM'):
                    ligand = pdb_line[17:26].strip()
                    this_ligand_charge = pdb_line[77:80].strip()
                    if ligand == ligand_name and this_ligand_charge in charge_dict:
                        ligand_charges[ligand_name] += int(charge_dict[this_ligand_charge])
        for ligand_name, charge in ligand_charges.items():
            for i, line in enumerate(log_lines):
                if line.startswith("Chain"):
                    continue
                elif ligand_name in line:
                    log_lines[i] = line.strip() + f' charge {charge}\n'
        with open("FMOPhore.log", 'w') as log_file:
            log_file.writelines(log_lines)
        return pdb_lines
    ################################################################
    def sep_chains(self, ligand_info, multiple_chains=True):
        if multiple_chains == False:
            output_folder = f"{self.pdb_file[:-4]}_{ligand_info}"
            os.makedirs(output_folder, exist_ok=True)
            io = PDBIO()
            pdb_name = os.path.basename(f"{self.pdb_file[:-4]}_{ligand_info}.pdb")
            pdb_id = os.path.splitext(pdb_name)[0]
            try:
                PDB_file = PDBParser().get_structure(pdb_id, f"{self.pdb_file[:-4]}_{ligand_info}.pdb")
                for chain in PDB_file.get_chains():
                    if f"{chain.id}" == ligand_info[4]:
                        ligand_file_path = os.path.join(output_folder, f"{pdb_id}.pdb")
                        io.set_structure(chain)
                        io.save(ligand_file_path)
            except Exception as e:
                print(f"{pdb_id}: {e}", file=open("FMOPhore.log", "a"))
        elif multiple_chains == True:
            output_folder = f"{self.pdb_file[:-4]}_{ligand_info}"
            os.makedirs(output_folder, exist_ok=True)
        new_pdb_path = os.path.join(output_folder, f"{self.pdb_file[:-4]}_{ligand_info}.pdb")
        shutil.move(f"{self.pdb_file[:-4]}_{ligand_info}.pdb", new_pdb_path) 
    ################################################################
    def rename_LIG(self, pdb_lines):
        excluded_strings = ['FMT', 'IOD', 'HOH', 'WAT', 'T3P', 'CL','NAG' , 'PO4',
                            'MG', 'ZN','PEG', 'DMS', 'GOL', 'BO3',  'NO3',
                            'EDO', 'SO4', 'ACE', 'NMA', 'NME', 'ACT', 'TRS',
                            'MES', 'OCY','FMT','PEG','SEP', 'BME', 'CSO', 
                            'IMD','TPO', 'TCE', 'MDP', 'NI'
                            , 'SAH', 'MTA', 'MG'
                            ]
        for i in range(len(pdb_lines)):
            if pdb_lines[i][17:20] == "LIG" and pdb_lines[i][0:6] != "TER   ":
                pdb_lines[i] = pdb_lines[i][:0] + "HETATM" + pdb_lines[i][6:]
        for i, line in enumerate(pdb_lines):
            if line.startswith('HETATM') and line[17:20].strip() not in excluded_strings:
                pdb_lines[i] = line[:22] + '9999' + line[26:]
        for i, line in enumerate(pdb_lines):
            record_type = line[:6].strip()
            residue_name = line[17:20].strip()
            if record_type in ['ATOM', 'HETATM']:
                if record_type == 'HETATM':
                    if line[17:20].strip() in ['T3P', 'WAT']:
                        line = line[:17] + 'HOH' + line[20:]
                    elif line[17:20].strip() in ['HIE', 'HIP', 'HID']:
                        line = line[:17] + 'HIS' + line[20:]
                if record_type == 'HETATM':
                    if line[17:20].strip() not in excluded_strings:
                        line = line[:17] + 'LIG' + line[20:]
                    if line[26].isdigit() or line[26].isalpha():
                        line = line[:26] + ' ' + line[26:]
                if record_type == 'ATOM':
                    if line[16:17].strip() in ['A']:
                        line = line[:16] + ' ' + line[17:]
                if record_type in ['ATOM', 'HETATM']:
                    if line[73:74].strip() in ['A', 'B', 'C', 'D']:
                        line = line[:73] + ' ' + line[74:]
                pdb_lines[i] = line
        return pdb_lines
    ################################################################
    def modify_atom_names(self, pdb_lines):
        atom_names = set()
        for i, line in enumerate(pdb_lines):
            record_type = line[:6].strip()
            if record_type in ['HETATM']:
                residue_name = line[17:20].strip()
                atom_name = line[11:16].strip()
                if residue_name + atom_name in atom_names:
                    if atom_name[-1].isdigit():
                        new_atom_name = atom_name[:-1] + str(int(atom_name[-1]) + 3)
                    else:
                        new_atom_name = atom_name
                    pdb_lines[i] = line[:12] + new_atom_name.ljust(4) + line[16:]
                else:
                    atom_names.add(residue_name + atom_name)
        return pdb_lines
    ################################################################
    def save_processed_pdb_file_ligand_info(self, pdb_lines, ligand_info):
        with open(f"{self.pdb_file[:-4]}_{ligand_info}.pdb", 'w') as f:
            f.writelines(pdb_lines)
    ################################################################
    def renumber_water_molecules(self, ligand_info):
        pdb_file_path = self.pdb_file[:-4] + f'_{ligand_info}.pdb'
        water_residues = []
        other_residues = []
        ligand_res_number = None
        with open(pdb_file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith(("ATOM", "HETATM")):
                    res_name = line[17:20].strip()
                    res_number = int(line[22:26].strip())
                    if res_name == 'LIG':
                        ligand_res_number = res_number
                    if res_name == 'HOH':
                        water_residues.append(line)
                    else:
                        other_residues.append(line)
        if ligand_res_number is None:
            raise ValueError("Ligand (LIG) not found in the PDB file.")
        new_hoh_number = ligand_res_number - 1
        renumbered_waters = []
        current_residue_number = None
        water_count = 0
        for line in water_residues:
            res_number = int(line[22:26].strip())  
            if current_residue_number is None or res_number != current_residue_number:
                water_count += 1
                new_res_id = new_hoh_number - water_count  
                current_residue_number = res_number  
            if new_res_id < 1:
                new_res_id = 1
            new_line = line[:22] + f"{new_res_id:4d}" + line[26:]
            renumbered_waters.append(new_line)
        output_filename = self.pdb_file[:-4] + f'_{ligand_info}.pdb'
        with open(output_filename, 'w') as out_file:
            for residue in other_residues:
                out_file.write(residue)
            for residue in renumbered_waters:
                out_file.write(residue)
    def renumber_important_residues(self, ligand_info):
        pdb_file_path = self.pdb_file[:-4] + f'_{ligand_info}.pdb'
        important_residues = ['HOH', 'ZN', 'CL', 'NI', 'CO', 'MG']
        renumbered_residues = []
        other_residues = []
        ligand_res_number = None
        with open(pdb_file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith(("ATOM", "HETATM")):
                    res_name = line[17:20].strip()
                    res_number = int(line[22:26].strip())
                    if res_name == 'LIG':
                        ligand_res_number = res_number
                    if res_name in important_residues:
                        renumbered_residues.append(line)
                    else:
                        other_residues.append(line)
        if ligand_res_number is None:
            raise ValueError("Ligand (LIG) not found in the PDB file.")
        new_res_number = ligand_res_number - 1
        renumbered_lines = []
        current_residue_number = None
        residue_count = 0
        for line in renumbered_residues:
            res_number = int(line[22:26].strip())
            if current_residue_number is None or res_number != current_residue_number:
                residue_count += 1
                new_res_id = new_res_number - residue_count
                current_residue_number = res_number
            if new_res_id < 1:
                new_res_id = 1
            new_line = line[:22] + f"{new_res_id:4d}" + line[26:]
            renumbered_lines.append(new_line)
        output_filename = self.pdb_file[:-4] + f'_{ligand_info}.pdb'
        with open(output_filename, 'w') as out_file:
            for residue in other_residues:
                out_file.write(residue)
            for residue in renumbered_lines:
                out_file.write(residue)
###################################################################################################
if __name__ == "__main__":
    """
    FMOPhore V.0.1 - LIGProcessor - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    LIG_processor = LIGProcessor(args.pdb_file)
    LIG_processor.process_complex_LIG()
