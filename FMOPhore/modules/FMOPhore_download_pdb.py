import argparse
import urllib.error
import os
from urllib.request import urlretrieve
from .FMOPhore_utility import EnvironmentGuard
# EnvironmentGuard().enforce()

"""
FMOPhore V.0.1 - PDB downloader - Copyright "©" 2024, Peter E.G.F. Ibrahim.
"""

def download_pdb(pdb_id, error_pdb_ids=None, dest_dir="."):
    """Download a single PDB file from RCSB into the specified directory."""
    if error_pdb_ids is None:
        error_pdb_ids = []
    pdb_id = pdb_id.strip()
    if not pdb_id:
        return
    pdb_file = f"{pdb_id}.pdb"
    pdb_url = f"https://files.rcsb.org/download/{pdb_file}"
    dest_path = os.path.join(dest_dir, pdb_file)
    try:
        urlretrieve(pdb_url, dest_path)
        print(f"Downloaded: {dest_path}")
    except urllib.error.HTTPError:
        print(f"Failed to download: {pdb_id}")
        error_pdb_ids.append(pdb_id)

def download_pdbs_from_file(file_path, error_pdb_ids=None, dest_dir="."):
    """Download multiple PDB files listed in a file into the specified directory."""
    if error_pdb_ids is None:
        error_pdb_ids = []
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            pdb_id = line.strip()
            if pdb_id:  # skip empty lines
                download_pdb(pdb_id, error_pdb_ids, dest_dir=dest_dir)
    if error_pdb_ids:
        print("\nThe following PDB IDs could not be downloaded:")
        for failed_id in error_pdb_ids:
            print(failed_id)
