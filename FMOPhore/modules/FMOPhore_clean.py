import numpy as np
import seaborn as sns
import pandas as pd
import os 
from itertools import islice
from math import factorial
import matplotlib.pyplot as plt
import re
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import csv
from numpy import nan
import glob
from rdkit import Chem
from rdkit.Chem import Draw
import argparse
from tabulate import tabulate
from PIL import Image, ImageDraw, ImageFont
import math
from .FMOPhore_utility import EnvironmentGuard
# EnvironmentGuard().enforce()

class Clean:
    def __init__(self, pdb_file, Binding_Energy=False, FE_score=False):
        self.pdb_file = pdb_file
        self.Binding_Energy = Binding_Energy
        self.FE_score = FE_score

    def Cleaning(self):
        ###################################################################################################
        if self.Binding_Energy:
            os.remove(f'../analysis/{self.pdb_file[:-4]}_Binding_Energy.txt')
            os.remove(f'../analysis/{self.pdb_file[:-4]}_complex_energy.txt')
            os.remove(f'../analysis/{self.pdb_file[:-4]}_protein_energy.txt')
            os.remove(f'../analysis/{self.pdb_file[:-4]}_ligand_energy.txt')
            os.remove(f'../analysis/{self.pdb_file[:-4]}_Affinity.txt')
        
        if self.FE_score:
            prefix = extract_prefix(self.pdb_file)
            print(f"Prefix extracted: {prefix}")
            file_patterns = [
                f'../analysis/{prefix}FMO_lig_H_frag*_H_FE_PIEDA_values.txt',
                f'../analysis/{prefix}FMO_lig_H_frag*_H_FE_FMO.txt',
                f'../analysis/{prefix}FMO_lig_H_frag*_H_FE_total_average.txt'
            ]
            for pattern in file_patterns:
                matched_files = glob.glob(pattern)
                for file in matched_files:
                    try:
                        os.remove(file)
                    except OSError as e:
                        print(f"Error removing {file}: {e}")

        common_files = [
            f'../analysis/{self.pdb_file[:-4]}_PIEDA_values.txt',
            f'../analysis/{self.pdb_file[:-4]}_FMO.txt',
            f'../analysis/{self.pdb_file[:-4]}_total_average.txt',
            f'Ph4_3D.txt'
        ]
        for file in common_files:
            if os.path.exists(file):
                try:
                    os.remove(file)
                except OSError as e:
                    print(f"Error removing {file}: {e}")



                        # try:
                        #     print(f"analysis -com ../{self.pdb_file[:-4]}_{ligand_info}_FMO_{distance}_H.pdb")
                        #     log_files = glob.glob('*.log')
                        #     inp_files = glob.glob('*.inp')
                        #     pdb_files = [f for f in glob.glob('*.pdb')]
                        #     files_to_clean = log_files + inp_files + pdb_files
                        #     for file in files_to_clean:
                        #         if os.path.exists(file):
                        #             try:
                        #                 os.remove(file)
                        #                 print(f"Removed: {file}")
                        #             except OSError as e:
                        #                 print(f"Error removing {file}: {e}")
                        #         else:
                        #             print(f"File not found: {file}")

                        #     os.chdir("..")
                        #     clean_up(['DFTB', 'rungms', 'rungms', 'rungms_personalized', 'DFTB_personalized'])
                        #     print("parent_dir = " , os.getcwd(), "aaaaaaaa")
                        #     os.chdir('..')
                        # except OSError as e:
                        #     if e.errno == 116:
                        #         print(f"Error changing directory: {e}")
                        #     else:
                        #         pass
                        # except Exception as e:
                        #     pass


if __name__ == "__main__":
    """
    FMOPhore V 0.1 - Analysis - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    parser = argparse.ArgumentParser(description="FMOPhore Fragmention Analysis")
    parser.add_argument("-pdb", "--pdb_file", type=str, required=True, help="Path to the PDB file")
    parser.add_argument('-FE', '--FE_score', action='store_true', default=None, help=argparse.SUPPRESS)
    parser.add_argument('-BE', '--Binding_Energy', action='store_true', default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()
    pdb_processor = Clean(args.pdb_file, args.Binding_Energy, args.FE_score)
    pdb_processor.Cleaning()
