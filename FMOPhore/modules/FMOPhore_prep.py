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
import argparse
import string
from .FMOPhore_chain_correction import chain_corrector
from collections import defaultdict, Counter
from .FMOPhore_utility import EnvironmentGuard
# EnvironmentGuard().enforce()

###################################################################################################
#  PDBProcessor
###################################################################################################
class PDBProcessor:
    def __init__(self, pdb_file, cofactor):
        self.pdb_file = pdb_file
        self.cofactor = cofactor
    def chain_id_list_all(self):
        return self.chain_id_list
    def ligands_list_all(self):
        return self.ligands_list
    ################################################################
    chain_id_list = []
    ligands_list  = []
    def process_complex(self, cofactor):
        pdb_lines = self.load_pdb_file()
        pdb_lines = self.filter_lines(pdb_lines)
        pdb_lines = self.renumber_residues_incrementally(pdb_lines)
        pdb_lines = self.residue_std(pdb_lines)
        pdb_lines = self.rename_LIG(pdb_lines)
        pdb_lines = self.increment_HETATM(pdb_lines)
        pdb_lines = self.modify_LIG_atom_names(pdb_lines)
        pdb_lines = self.recenter(pdb_lines)
        pdb_lines = self.load_pdb_file()
        pdb_lines = self.chain_id(pdb_lines)
        pdb_lines = self.write_cofactor(pdb_lines, cofactor)
        # self.save_pdb_file(pdb_lines)
        pdb_lines = self.load_pdb_file()
        pdb_lines = self.rechain_LIG(pdb_lines)
        pdb_lines = self.set_ligand_chain(pdb_lines)
        pdb_lines = self.correct_ligand_chain(pdb_lines)
        pdb_lines = self.remove_duplicate_water_molecules(pdb_lines)
        pdb_lines = self.ligand_charge(pdb_lines)
        # # pdb_lines = self.renumber_atoms(pdb_lines)
        self.save_pdb_file(pdb_lines)
        pdb_lines = self.pdbfixer()
    ################################################################
    def load_pdb_file(self):
        with open(self.pdb_file, 'r') as f:
            return f.readlines()
    ###############################################################
    def filter_lines(self, pdb_lines):
        excluded_strings = sorted(set([
            'HD2 ASP', 'HA2 PHE', 'FMT', 'ANISOU', 'BU3', 'NAG', 'PO4',
            'PEG', 'GOL', 'BO3', 'EDO', 'SO4', 'NO3', 'ACE', 'NMA', 'NME',
            'ACT', 'MES', 'TRS', 'OCY', 'SEP', 'BME', 'CSO', 'IMD', 'TPO',
            'TCE', 'MDP', 'IOD', 'IPA', 'PTR',
            'T3P', 'WAT', 'DMS', 'EOH', 'GZ6', 'PCA', 'ACY', 'MLA'
        ]))
        excluded_metals = {'CA', 'NI', 'CO', 'ZN', 'MG'}
        filtered_lines = []
        for line in pdb_lines:
            record_type = line[0:6].strip()
            resname = line[17:20].strip()
            chain_id = line[21:22]
            altloc = line[16:17]
            atom_name = line[12:16].strip()
            if (
                any(substr in line for substr in excluded_strings) or
                resname in ['CD', 'CL', 'NA'] or
                altloc in ['B', 'C', 'D', 'E', 'F'] or
                # line[26:27] == 'A' or
                (record_type == "ATOM" and atom_name == 'OXT') or
                (record_type == "HETATM" and resname in excluded_metals)
            ):
                continue
            if line[16:17] == 'A':
                line = line[:16] + ' ' + line[17:]
            # All uppercase alphabetic insertion codes
            insertion_codes_to_strip = set(string.ascii_uppercase)
            if line[26:27] in insertion_codes_to_strip:
                line = line[:26] + ' ' + line[27:]
            filtered_lines.append(line)

        for i, line in enumerate(filtered_lines):
            record_type = line[:6].strip()
            residue_name = line[17:20].strip()
            if record_type == 'HETATM' and residue_name != 'HOH':
                if line[22].isdigit() or line[22].isalpha():
                    line = line[:22] + ' ' + line[23:]
            if line.startswith('HETATM') and residue_name not in excluded_strings and residue_name != 'HOH':
                residue_name = line[17:20].strip()
                resnum_field = line[22:26]
                resnum_clean = resnum_field.strip()
                # Replace numeric residue name with 'LIG'
                if residue_name.isdigit():
                    line = line[:17] + 'LIG' + line[20:]
                # Replace 1-digit or 4-digit residue number with 999
                if resnum_clean.isdigit() and (len(resnum_clean) == 1 or len(resnum_clean) == 4):
                    line = line[:22] + '{:>4}'.format('999') + line[26:]
                filtered_lines[i] = line  # Update line after both modifications       
        return filtered_lines
    ###############################################################
    def residue_std(self, pdb_lines):
        residue_std_lines = []
        current_residue = None
        hz_found = set()
        for i, line in enumerate(pdb_lines):
            if line[17:20].strip() == 'LYS':  
                residue_number = line[22:26].strip()
                atom_name = line[12:16].strip()
                if residue_number != current_residue:
                    current_residue = residue_number
                    hz_found.clear()
                if atom_name in ['HZ1', 'HZ2', 'HZ3']:
                    hz_found.add(atom_name)
                if atom_name == 'HZ':
                    if 'HZ1' not in hz_found:
                        new_atom_name = 'HZ1'
                    elif ' HZ2' not in hz_found:
                        new_atom_name = 'HZ2'
                    else:
                        new_atom_name = 'HZ3'
                    pdb_lines[i] = line[:12] + f"{new_atom_name:<4}" + line[16:]  # Replace atom name in line
        return pdb_lines
    ###############################################################
    def renumber_residues_incrementally(self, pdb_lines):
        corrected_lines = []
        prev_res_num = None
        prev_res_name = None
        current_chain_id = None
        res_offset = 0 
        for line in pdb_lines:
            if not line.startswith(("ATOM", "HETATM")):
                corrected_lines.append(line)
                continue
            res_name = line[17:20].strip()
            chain_id = line[21]
            res_num_str = line[22:26]
            res_num = int(res_num_str.strip())
            if chain_id != current_chain_id:
                current_chain_id = chain_id
                prev_res_num = res_num
                prev_res_name = res_name
                res_offset = 0
                new_res_num = res_num
            elif res_num == prev_res_num and res_name != prev_res_name:
                res_offset += 1
                new_res_num = res_num + res_offset
            else:
                new_res_num = res_num + res_offset
            prev_res_num = res_num
            prev_res_name = res_name
            new_res_str = f"{new_res_num:>4}"
            new_line = line[:22] + new_res_str + line[26:]
            corrected_lines.append(new_line)
        return corrected_lines

    ################################################################
    def rename_LIG(self, pdb_lines):
        def is_float(value):
            try:
                float(value)
                return True
            except ValueError:
                return False
        excluded_mistaken_errors = [
                                    "TER", "HEADER", "REMARK", 
                                    "END", "CRYST1", "JRNL", 
                                    "SOURCE", "COMPND", 
                                    "KEYWDS", "SSBOND",
                                    "REVDAT", "FORMUL",
                                    "SITE", "CISPEP",
                                    "DBREF","SEQRES", "HET", 
                                    "HELIX", "SHEET",
                                    "SSBOND", "CISPEP", 
                                    "ORIGX1", "SCALE1", 
                                    "ORIGX2", "SCALE2",
                                    "ORIGX3", "SCALE3",  
                                    "MTRIX1", "LINK",
                                    "MTRIX2", "HETSYN",
                                    "MTRIX3", "AUTHOR",
                                    "TITLE", "MASTER",
                                    "HETNAM", "FORMUL", 
                                    "CONECT", "SEQADV",
                                    "EXPDTA", "MODRES",
                                    "CRYST1", "ANISOU", "MODEL"
                                    ]
        residues_names = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
                          'GLN', 'GLU', 'GLY', 'HIS', 'HIE', 
                          'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                          'PRO', 'SER', 'THR', 'TRP', 'TYR', 
                          'VAL']
        for i in range(len(pdb_lines)):
            if (pdb_lines[i][17:20].strip() not in residues_names and 
                pdb_lines[i][6:12].strip().isdigit() and 
                is_float(pdb_lines[i][30:38].strip()) and 
                pdb_lines[i][:6].strip() not in excluded_mistaken_errors):  
                pdb_lines[i] = pdb_lines[i][:0] + "HETATM" + pdb_lines[i][6:]
        return pdb_lines
    ################################################################
    def increment_HETATM(self, pdb_lines):
        last_chain = None
        base_increment = 0
        output_lines = []
        for line in pdb_lines:
            if line.startswith("HETATM"):  
                chain = line[21]
                residue = int(line[22:26].strip())
                if chain != last_chain:
                    if last_chain is not None:
                        base_increment += 0  
                    last_chain = chain
                new_residue = residue + base_increment
                output_line = line[:22] + f"{new_residue:4}" + line[26:]
                output_lines.append(output_line)
            else:
                output_lines.append(line)
        return output_lines
    ###############################################################    
    def modify_LIG_atom_names(self, pdb_lines):
        atom_names = set()
        for i, line in enumerate(pdb_lines):
            record_type = line[:6].strip()
            if record_type in ['HETATM']:
                residue_number = line[22:26].strip()
                residue_name = line[17:26].strip()
                atom_name = line[12:16].strip()
                unique_identifier = residue_name + atom_name
                if unique_identifier in atom_names:
                    if atom_name[-1].isdigit():
                        increment = 1
                        while unique_identifier in atom_names:
                            num = ''.join(filter(str.isdigit, atom_name))
                            new_num = str(int(num) + increment)
                            new_atom_name = atom_name.rstrip('0123456789') + new_num
                            new_atom_name = new_atom_name[:4] 
                            unique_identifier = residue_name + new_atom_name
                            increment += 1
                    else:
                        new_atom_name = atom_name[:3] + '1'  
                        unique_identifier = residue_name + new_atom_name
                    pdb_lines[i] = line[:12] + new_atom_name.ljust(4) + line[16:]
                atom_names.add(unique_identifier)
        return pdb_lines
    ################################################################
    def recenter(self, pdb_lines):
        atom_coords = []
        atom_masses = []
        for line in pdb_lines:
            if line.startswith(("ATOM", "HETATM")): 
                atom_type = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                mass = float(line[54:60]) if line[54:60].strip() else 1.0
                atom_coords.append([x, y, z])
                atom_masses.append(mass)
        atom_coords = np.array(atom_coords)
        atom_masses = np.array(atom_masses)
        total_mass = np.sum(atom_masses)
        center_of_mass = np.sum(atom_coords * atom_masses[:, np.newaxis], axis=0) / total_mass
        recentered_coords = atom_coords - center_of_mass
        new_pdb_lines = []
        line_index = 0
        for line in pdb_lines:
            if line.startswith(("ATOM", "HETATM")):
                atom_type = line[12:16].strip()
                x, y, z = recentered_coords[line_index]
                new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
                if not new_line.endswith('\n'):
                    new_line += '\n'
                new_pdb_lines.append(new_line)
                line_index += 1
        output_file = self.pdb_file
        with open(output_file, 'w') as f:
            for line in new_pdb_lines:
                if 'TER' not in line and line.strip():
                    f.write(line)
        return new_pdb_lines
    ################################################################        
    def chain_id(self, pdb_lines):
        chain_ranges = {}
        current_chain_id = None
        current_range_start = None
        current_range_end = None
        for line in pdb_lines:
            if line.startswith("ATOM"):
                chain_id = line[21]
                residue_number = int(line[22:26])
                if chain_id != current_chain_id:
                    if current_chain_id != " " and current_chain_id is not None:
                        chain_ranges[current_chain_id] = (current_range_start, current_range_end)
                    current_chain_id = chain_id
                    current_range_start = residue_number
                current_range_end = residue_number
        if current_chain_id != " ":
            chain_ranges[current_chain_id] = (current_range_start, current_range_end)
            with open("FMOPhore.log", "w") as file:
                for chain_id, (start, end) in chain_ranges.items():
                    file.write(f"Chain {chain_id}: {start}-{end}\n")
        if current_chain_id == " ":
            chain_ranges = {}
            with open("../FMOPhore_org.log", "r") as file:
                lines = file.readlines()
                with open("FMOPhore.log", "w") as new_file:
                    lines_with_newline = [line.rstrip("\n") + "\n" for line in lines]
                    new_file.writelines(lines_with_newline)
                for line in lines:
                    if line.startswith("Chain"):
                        chain_id, range_str = line.strip().split(": ")
                        start_residue, end_residue = map(int, range_str.split("-"))
                        chain_ranges[chain_id] = (start_residue, end_residue)

            chain_ranges_str = [f"{start-1}-{end}" for _, (start, end) in chain_ranges.items()]
            error_message = "Running chain correction!"
            with open("errors.log", "a") as error_file:
                error_file.write(error_message + "\n")
            chain_corrector(pdb_file_pattern=self.pdb_file, chain_from_to=chain_ranges_str)
    ################################################################        
    def write_cofactor(self, pdb_lines, cof: dict):
        if pdb_lines is None:
            with open(self.pdb_file, "r", encoding="utf-8", errors="ignore") as fh:
                pdb_lines = fh.readlines()
        if not cof or not cof.get("name"):
            return pdb_lines
        target_name = str(cof["name"]).upper().strip()
        target_resnum = None
        if cof.get("resnum") not in (None, ""):
            try:
                target_resnum = int(cof["resnum"])
            except (TypeError, ValueError):
                target_resnum = None
        candidates = []
        for line in pdb_lines:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:20].strip().upper()
            if resname != target_name:
                continue
            chain = line[21]
            resnum_str = line[22:26].strip()  
            if not resnum_str.isdigit():
                continue
            resnum = int(resnum_str)
            if target_resnum is None or resnum == target_resnum:
                candidates.append((resname, chain, resnum))
        if candidates:
            if target_resnum is not None:
                matches = [c for c in candidates if c[2] == target_resnum]
                chosen = matches[0] if matches else candidates[0]
            else:
                from collections import Counter
                (chain, resnum), _ = Counter((c[1], c[2]) for c in candidates).most_common(1)[0]
                chosen = (target_name, chain, resnum)
            name, chain, resnum = chosen
            chain_out = chain if chain != " " else "_"  
            with open("FMOPhore.log", "a", encoding="utf-8") as f:
                f.write(f"cofactor: {name} {chain_out} {resnum} !\n")
        else:
            with open("FMOPhore.log", "a", encoding="utf-8") as f:
                f.write(f"cofactor: {target_name}-_-NA  # not found in PDB\n")
        return pdb_lines
    ################################################################
    def rechain_LIG(self, pdb_lines):
        excluded_strings = sorted(set([
            'HD2 ASP', 'HA2 PHE', 'FMT', 'ANISOU', 'BU3', 'NAG', 'PO4',
            'PEG', 'GOL', 'BO3', 'EDO', 'SO4', 'NO3', 'ACE', 'NMA', 'NME',
            'ACT', 'MES', 'TRS', 'OCY', 'SEP', 'BME', 'CSO', 'IMD', 'TPO',
            'TCE', 'MDP', 'IOD', 'IPA', 'PTR',
            'T3P', 'WAT', 'MG', 'DMS', 'EOH', 'GZ6', 'PCA', 'ACY', 'MLA'
        ]))
        for i, hetatm_line in enumerate(pdb_lines):
            if hetatm_line.startswith('HETATM') and hetatm_line[17:20].strip() not in excluded_strings:
                hetatm_chain_id = hetatm_line[21]
                hetatm_residue_number = int(hetatm_line[22:26])
                hetatm_coords = [float(hetatm_line[30:38]), float(hetatm_line[38:46]), float(hetatm_line[46:54])]
                nearest_chain_id = None
                min_distance = float('inf')
                for atom_line in pdb_lines:
                    if atom_line.startswith("ATOM"):
                        atom_chain_id = atom_line[21]
                        atom_residue_number = int(atom_line[22:26])
                        atom_coords = [float(atom_line[30:38]), float(atom_line[38:46]), float(atom_line[46:54])]
                        distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(hetatm_coords, atom_coords)]))
                        if distance < min_distance:
                            min_distance = distance
                            nearest_chain_id = atom_chain_id
                if nearest_chain_id:
                    pdb_lines[i] = hetatm_line[:21] + nearest_chain_id + hetatm_line[22:]         
        return pdb_lines
    # ################################################################
    def set_ligand_chain(self, pdb_lines):
        excluded_strings = sorted(set([
            'HD2 ASP', 'HA2 PHE', 'FMT', 'ANISOU', 'BU3', 'NAG', 'PO4',
            'PEG', 'GOL', 'BO3', 'EDO', 'SO4', 'NO3', 'ACE', 'NMA', 'NME',
            'ACT', 'MES', 'TRS', 'OCY', 'SEP', 'BME', 'CSO', 'IMD', 'TPO',
            'TCE', 'MDP', 'IOD', 'IPA', 'PTR',
            'T3P', 'WAT', 'MG', 'DMS', 'EOH', 'GZ6', 'PCA', 'ACY', 'MLA'
        ]))
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
                            if ligand_key[0].upper() not in [' CL', 'T3P', 'WAT', 'HOH']:
                                error_message = "{} falls between multiple chains {} !".format(ligand, chain_id)
                                error_file.write(error_message + "\n")
                        else:
                            if ligand_key[0].upper() not in [' CL', 'T3P', 'WAT', 'HOH']:
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
            if record_type == 'HETATM' and line[17:20].strip() not in [' CL','T3P', 'WAT', 'HOH']:
                chain_id = line[21:22]
                ligand = line[17:26]
                self.chain_id_list.append(chain_id)  
                self.ligands_list.append(ligand)
        return pdb_lines   
    ################################################################
    def correct_ligand_chain(self, pdb_lines):
        excluded_strings = sorted(set([
            'HD2 ASP', 'HA2 PHE', 'FMT', 'ANISOU', 'BU3', 'NAG', 'PO4',
            'PEG', 'GOL', 'BO3', 'EDO', 'SO4', 'NO3', 'ACE', 'NMA', 'NME',
            'ACT', 'MES', 'TRS', 'OCY', 'SEP', 'BME', 'CSO', 'IMD', 'TPO',
            'TCE', 'MDP', 'IOD', 'IPA', 'PTR',
            'T3P', 'WAT', 'MG', 'DMS', 'EOH', 'GZ6', 'PCA', 'ACY', 'MLA'
        ]))

        self.ligands_list.clear()
        self.chain_id_list.clear()

        ligand_chain_map = defaultdict(list)

        with open("FMOPhore.log", "r") as log_file:
            for line in log_file:
                line = line.strip()
                if "falls in chain " in line:
                    parts = line.split()
                    if len(parts) >= 5:
                        ligand_name = parts[0].strip()
                        chain_id = parts[1].strip()
                        ligand_num = parts[2].strip()
                        ligand_key = (ligand_name, ligand_num)
                        ligand_chain_map[ligand_key].append(chain_id)
        ligand_to_chain = {
            key: Counter(chains).most_common(1)[0][0]
            for key, chains in ligand_chain_map.items()
        }
        seen = set()
        for i, line in enumerate(pdb_lines):
            ligand_name = line[17:20].strip()
            if line.startswith("HETATM") and ligand_name not in excluded_strings:
                ligand_name = line[17:20].strip()
                ligand_num = line[22:26].strip()
                chain_id = line[21]
                ligand_key = (ligand_name, ligand_num)
                if ligand_key in ligand_to_chain:
                    correct_chain = ligand_to_chain[ligand_key]
                    if chain_id != correct_chain:
                        pdb_lines[i] = line[:21] + correct_chain + line[22:]
                        chain_id = correct_chain  
                ligand_id = f"{ligand_name} {chain_id} {ligand_num}"
                if ligand_id not in seen:
                    seen.add(ligand_id)
                    self.chain_id_list.append(chain_id)
                    self.ligands_list.append(ligand_id)
        return pdb_lines

    ################################################################
    def remove_duplicate_water_molecules(self, pdb_lines):
        clean_pdb_lines = []
        seen_water_residues = set()
        
        for line in pdb_lines:
            if line.startswith('HETATM') and line[17:20].strip() == 'HOH':
                chain_id = line[21].strip()
                residue_number = line[22:26].strip()
                atom_name = line[12:16].strip()
                atom_cords = line[0:55].strip()
                residue_key = (chain_id, residue_number, atom_name, atom_cords)

                if residue_key not in seen_water_residues:
                    seen_water_residues.add(residue_key)
                    clean_pdb_lines.append(line)
                else:
                    # Skip duplicate HOH residue atoms
                    continue
            elif line.startswith('TER'):
                # Keep TER only if it's for a water that was preserved
                if (
                    len(clean_pdb_lines) > 0 and
                    clean_pdb_lines[-1].startswith('HETATM') and
                    clean_pdb_lines[-1][17:20].strip() == 'HOH'
                ):
                    clean_pdb_lines.append(line)
            else:
                clean_pdb_lines.append(line)
        
        return clean_pdb_lines
    ################################################################
    def ligand_charge(self, pdb_lines):
        excluded_strings = sorted(set([
            'HD2 ASP', 'HA2 PHE', 'FMT', 'ANISOU', 'BU3', 'NAG', 'PO4',
            'PEG', 'GOL', 'BO3', 'EDO', 'SO4', 'NO3', 'ACE', 'NMA', 'NME',
            'ACT', 'MES', 'TRS', 'OCY', 'SEP', 'BME', 'CSO', 'IMD', 'TPO',
            'TCE', 'MDP', 'IOD', 'NI', 'IPA', 'CO', 'ZN', 'PTR',
            'T3P', 'MG', 'DMS', 'EOH', 'GZ6', 'PCA', 'ACY', 'MLA'
        ]))
        charge_dict = {'O1+': +1, 'N1+': +1, 'O1-': -1, 'S1-': -1, 'l1-': -1, 'N1-': -1}
        ligand_charges = {}
        for line in pdb_lines:
            if line.startswith("HETATM") and line[17:20].strip() not in excluded_strings:
                resname = line[17:20].strip()
                chain = line[21].strip()
                resnum = line[22:26].strip()
                key = f"{resname} {chain} {resnum}"
                charge_label = line[76:80].strip()
                if charge_label in charge_dict:
                    ligand_charges.setdefault(key, 0)
                    ligand_charges[key] += charge_dict[charge_label]
        with open("FMOPhore.log", 'r') as f:
            log_lines = f.readlines()
        updated_lines = []
        for line in log_lines:
            stripped = line.strip()
            if "Chain" in line or "multiple" in line:
                updated_lines.append(line)
                continue
            ligand_key = None
            for key in ligand_charges:
                if key in line:
                    ligand_key = key
                    break
            if "charge" not in line:
                if ligand_key:
                    updated_lines.append(stripped + f" charge {ligand_charges[ligand_key]}\n")
                else:
                    updated_lines.append(stripped + " charge 0\n")
            else:
                updated_lines.append(line)
        with open("FMOPhore.log", 'w') as f:
            f.writelines(updated_lines)
        return pdb_lines
    ################################################################
    def renumber_atoms(self, pdb_lines):
        updated_lines = []
        atom_counters = {}
        for line in pdb_lines:
            # if line.startswith("ATOM") or line.startswith("HETATM"):
            if line.startswith("HETATM"):
                residue_name = line[17:20].strip()
                residue_number = line[22:26].strip()
                chain_id = line[21].strip()
                atom_name = line[12:16].strip()
                residue_key = (residue_name, residue_number, chain_id)
                if residue_key not in atom_counters:
                    atom_counters[residue_key] = {}
                if atom_name.startswith("HX"):
                    base_name = "H"
                if atom_name[0].isdigit():
                    base_name = ''.join([char for char in atom_name[1:] if not char.isdigit()])
                else:
                    base_name = ''.join([char for char in atom_name if not char.isdigit()])
                if base_name not in atom_counters[residue_key]:
                    atom_counters[residue_key][base_name] = 0
                atom_counters[residue_key][base_name] += 1
                new_atom_name = f"{base_name}{atom_counters[residue_key][base_name]}"
                new_atom_name = new_atom_name[:4].ljust(4)
                updated_line = line[:12] + new_atom_name + line[16:]
                updated_lines.append(updated_line)
            else:
                updated_lines.append(line)
        return updated_lines
    ################################################################
    def save_pdb_file(self, pdb_lines):
        with open(f"{self.pdb_file[:-4]}.pdb", 'w') as f:
            f.writelines(pdb_lines)
    ################################################################
    def pdbfixer(self):
        original_pdb_path = f"{self.pdb_file[:-4]}.pdb"
        temp_fixed_path = f"{self.pdb_file[:-4]}_fixed.pdb"
        with open(original_pdb_path, "r") as file:
            original_lines = file.readlines()
        hoh_lines = [line for line in original_lines if line.startswith("HETATM") and line[17:20] == "HOH"]
        hoh_ter_lines = [line for line in original_lines if line.startswith("TER") and "HOH" in line]
        hetatm_non_hoh = [line for line in original_lines if line.startswith("HETATM") and line[17:20] != "HOH"]
        conect_lines = [line for line in original_lines if line.startswith("CONECT")]
        ter_lines = [line for line in original_lines if line.startswith("TER") and "HOH" not in line]
        end_line = [line for line in original_lines if line.strip() == "END"]
        os.system(
            f'pdbfixer {original_pdb_path} --output {temp_fixed_path} '
            f'--add-atoms=heavy --keep-heterogens=all --add-residues --ph=7.0'
        )
        with open(temp_fixed_path, "r") as file:
            fixed_lines = file.readlines()
        cleaned_lines = [line for line in fixed_lines if line.startswith("ATOM") or line.startswith("TER")]
        cleaned_lines.extend(hoh_lines)
        cleaned_lines.extend(hoh_ter_lines)
        cleaned_lines.extend(hetatm_non_hoh)
        amino_acid_mapping = {
            " NZ  LYS": "N1+",
            " OD2 ASP": "O1-",
            " NH2 ARG": "N1+",
            " OE2 GLU": "O1-"
        }
        his_line_count = {}
        for line in cleaned_lines:
            if "HIS" in line[17:20]:
                chain_id = line[21]
                residue_id = line[22:26]
                key = f"{chain_id}_{residue_id}"
                his_line_count[key] = his_line_count.get(key, 0) + 1
        for key, count in his_line_count.items():
            if count == 18:
                chain_id, residue_id = key.split("_")
                amino_acid_mapping[f" NE2 HIS{chain_id}{residue_id}"] = "N1+"
        final_lines = []
        for line in cleaned_lines:
            atom_key_short = line[12:20]
            atom_key_long = line[12:26]
            if atom_key_short in amino_acid_mapping:
                modified_line = line[:77] + amino_acid_mapping[atom_key_short] + line[80:]
                final_lines.append(modified_line)
            elif atom_key_long in amino_acid_mapping:
                modified_line = line[:77] + amino_acid_mapping[atom_key_long] + line[80:]
                final_lines.append(modified_line)
            else:
                final_lines.append(line)
        for line in ter_lines:
            if line not in final_lines:
                final_lines.append(line)
        final_lines.extend(conect_lines)
        if end_line:
            final_lines.append(end_line[0])
        else:
            final_lines.append("END\n")
        with open(original_pdb_path, "w") as file:
            file.writelines(final_lines)
        os.remove(temp_fixed_path)
###################################################################################################
if __name__ == "__main__":
    """
    FMOPhore V 0.1 - PDBProcessor - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    parser = argparse.ArgumentParser(description="FMOPhore Fragmention PDBProcessor")
    parser.add_argument("-pdb", "--pdb_file", type=str, required=True, help="Path to the PDB file")
    args = parser.parse_args()
    pdb_processor = PDBProcessor(args.pdb_file, args.cofactor)
    pdb_processor.process_complex(args.cofactor)