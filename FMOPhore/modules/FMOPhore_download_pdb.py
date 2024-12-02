import argparse
import urllib.error
from urllib.request import urlretrieve

"""
FMOPhore V.0.1 - PDB downloader - Copyright "©" 2024, Peter E.G.F. Ibrahim.
"""
def download_pdb(pdb_id, error_pdb_ids=None):
    if error_pdb_ids is None:
        error_pdb_ids = []
    pdb_file = f'{pdb_id}.pdb'
    pdb_url = f'https://files.rcsb.org/download/{pdb_file}'
    try:
        urlretrieve(pdb_url, pdb_file)
    except urllib.error.HTTPError as e:
        error_pdb_ids.append(pdb_id)
def download_pdbs_from_file(file_path, error_pdb_ids=None):
    if error_pdb_ids is None:
        error_pdb_ids = []
    pdb_ids = []
    with open(file_path, "r") as file:
        lines = file.readlines()
    for line in lines:
        pdb_id = line.strip().split()[0]
        pdb_ids.append(pdb_id)
    for pdb_id in pdb_ids:
        download_pdb(pdb_id, error_pdb_ids)