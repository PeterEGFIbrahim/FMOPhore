import os
from Bio.PDB import PDBParser
import numpy as np
import multiprocessing
import argparse

###################################################################################################
#  This part to define the distance and co-factors if any.
###################################################################################################
def distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2) 
###################################################################################################
#  Classes
###################################################################################################
excluded_residues = ['EOH', 'GZ6', 'CL','DMS', 'FMT', 
             'PEG', 'GOL', 'BO3', 'EDO', 
             'SO4', 'ACE', 'NMA', 'NME', 'ACT', 'MES', 
             'OCY', 'SEP', 'TPO', 'IPA', 'TRS', ' ZN', 'ZN']
class CutoffProcessor:
    def __init__(self, pdb_file, distance_cutoff, same_target=False):
        self.pdb_file = pdb_file
        self.distance_cutoff = distance_cutoff if distance_cutoff != "no_cutoff" else np.inf
        self.same_target = same_target
        self.chain_id_list = []
        self.ligands_list = []
    ################################################################
    def process_cutoff(self, ligand_info):
        if self.same_target:  # Correctly use self.same_target
            self.generate_cutoff_pdb_complex(distance_cutoff=self.distance_cutoff, ligand_info=ligand_info, same_target=True)
        else:
            self.generate_cutoff_pdb_complex(distance_cutoff=self.distance_cutoff, ligand_info=ligand_info)
        os.chdir("..")
    def cutoff_complex(self):
        pdb_lines = self.load_pdb_file()
        ligands_list = [self.pdb_file[-13:-4]]
        for ligand in ligands_list:
            self.process_cutoff(ligand)
    ################################################################
    def load_pdb_file(self):
        with open(self.pdb_file, 'r') as f:
            return f.readlines()
    ################################################################
    def generate_cutoff_pdb_complex(self, distance_cutoff, ligand_info, same_target=None):
        residues_to_select = []
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
                    # if residue.get_resname() not in excluded_residues:
                        ligand = residue
                        ligand_chain = chain.id
                        ligand_number = residue.id[1]
        residues = []
        selected_residues = set() 
        for model in renumbered_structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == ' ' or residue.get_id()[0] == 'H' or residue.get_resname() == 'HOH':
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

        os.remove(f'{self.pdb_file[:-4]}_residues.txt')        
        if same_target is True:
            with open(f'../../library_analysis/residues/{self.pdb_file[:-4]}_residues.txt', 'w') as file:
                file.write('\n'.join(updated_lines))
        else:
            with open(f'{self.pdb_file[:-4]}_residues.txt', 'w') as file:
                file.write('\n'.join(updated_lines))
################################################################################################################################
if __name__ == "__main__":
    """
    FMOPhore V 0.1 - CutoffProcessor - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    pdb_processor = CutoffProcessor(self.pdb_file)
    pdb_processor.cutoff_complex()
