import argparse
import glob
import os
import shutil

def renumber_pdb(input_filename, output_filename, chain_start, chain_end):
    current_chain_id = None
    current_residue_number = 0
    new_residue_number = 0

    with open(input_filename, 'r') as input_file:
        with open(output_filename, 'w') as output_file:
            for line in input_file:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain_id = line[21]
                    residue_number = int(line[22:26])
                    if chain_id != current_chain_id:
                        current_chain_id = chain_id
                        current_residue_number = 0
                    if residue_number != current_residue_number:
                        current_residue_number = residue_number
                        new_residue_number = 0
                    if residue_number <= chain_start:
                        new_residue_number = residue_number
                    elif chain_start < residue_number <= chain_end:
                        new_residue_number = residue_number - chain_start
                    else:
                        new_residue_number = residue_number 
                    modified_line = line[:22] + f'{new_residue_number:4d}' + line[26:]
                    output_file.write(modified_line)
                else:
                    output_file.write(line)

def assign_chain_id(chain_id, residue_number, previous_residue_number):
    if residue_number < previous_residue_number:
        chain_id = chr(ord(chain_id) + 1)
    return chain_id

def correct_chains(input_file, output_file):
    with open(input_file, 'r') as input_file:
        pdb_lines = input_file.readlines()
    current_chain_id = 'A'
    previous_residue_number = 0
    output_pdb_lines = []
    for line in pdb_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain_id = line[21]
            residue_number = int(line[22:26])
            current_chain_id = assign_chain_id(current_chain_id, residue_number, previous_residue_number)
            line = line[:21] + current_chain_id + line[22:]
            previous_residue_number = residue_number 
        output_pdb_lines.append(line)
    with open(output_file, 'w') as output_file_corrected:
        output_file_corrected.writelines(output_pdb_lines)
###################################################################################################
def chain_corrector(pdb_file_pattern=None, pdb_files_pattern=None, chain_from_to=None):
    """
    FMOPhore V 0.1 - Chain corrector - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    if not pdb_file_pattern and not pdb_files_pattern:
        raise ValueError("Provide either a PDB file path or a directory path.")
    if not chain_from_to:
        raise ValueError("Provide chain ranges.")
    if pdb_file_pattern:
        if os.path.isfile(pdb_file_pattern) and pdb_file_pattern.endswith(".pdb"):
            pdb_files = [pdb_file_pattern]
        else:
            raise ValueError("Invalid PDB file path or extension.")
    elif pdb_files_pattern:
        if os.path.isdir(pdb_files_pattern):
            pdb_files = glob.glob(os.path.join(pdb_files_pattern, "*.pdb"))
        else:
            raise ValueError("Invalid directory path.")
    else:
        raise ValueError("Provide valid file or directory inputs.")

    for input_path in pdb_files:
        for chain_range in chain_from_to:
            chain_start, chain_end = [int(x) for x in chain_range.split('-')]
            temp_output_filename = input_path + "_temp"
            renumber_pdb(input_path, temp_output_filename, chain_start, chain_end)
            shutil.move(temp_output_filename, input_path)
            correct_chains(input_path, input_path)
