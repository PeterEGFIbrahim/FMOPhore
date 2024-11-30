import os
import argparse
import glob

def merge_ligs_prot(protein_path, ligands_path, output_dir):
    """
    FMOPhore V.0.1 - Merger - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(protein_path, "rb") as protein_file:
        protein_lines = protein_file.readlines()
    protein_lines = [line for line in protein_lines if not (line.startswith(b"END") or line.startswith(b"CONECT"))]
    protein_contents = b"".join(protein_lines)
    with open(ligands_path, "rb") as ligand_file:
        ligand_lines = ligand_file.readlines()
    for index, line in enumerate(ligand_lines):
        if line.startswith(b"HETATM"):
            break
    ligand_contents = b"".join(ligand_lines[index:])
    ligand_filename = os.path.basename(ligands_path)
    output_path = os.path.join(output_dir, ligand_filename)
    with open(output_path, "wb") as output_file:
        output_file.write(protein_contents)
        output_file.write(ligand_contents)
