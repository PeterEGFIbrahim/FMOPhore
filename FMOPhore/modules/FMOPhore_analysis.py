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

class Analysis:
    def __init__(self, pdb_file, Binding_Energy=False, FE_score=False):
        self.pdb_file = pdb_file
        self.Binding_Energy = Binding_Energy
        self.FE_score = FE_score
    ###################################################################################################
    # create Ph4 + FMO file
    ###################################################################################################
    def analyze(self):
        def compare_files(file1, file2, output_file):
            with open(file1, 'r') as f1, open(file2, 'r') as f2, open(output_file, 'w') as output:
                lines1 = f1.readlines()
                lines2 = f2.readlines()
                file2_dict = {}
                for line2 in lines2[1:]:  
                    values2 = line2.split()
                    key = values2[0]
                    file2_dict[key] = line2.strip()
                output.write("Fragment\ttotal\tEes\tEex\tEct+mix\tEdisp\tGsol\tdistance\tAngle\tLIG.x\tLIG.y\tLIG.z\tPROT.x\tPROT.y\tPROT.z\tBond_type\tprot_atom_id\tlig_atom_id\n")
                # for line1 in sorted(lines1[1:], key=lambda x: int(x.split()[0].split('-')[1])):
                for line1 in lines1[1:]:
                    values1 = line1.split()
                    key = values1[0]
                    if key in file2_dict:
                        line2 = file2_dict[key]
                        formatted_values1 = [f"{float(val):8.3f}" for val in values1[1:-3]]
                        bond_type = values1[-3:-1]
                        output.write(line2 + '\t' + ''.join(formatted_values1) + '\t' + '\t'.join(bond_type) + '\t' + ''.join(values1[-1]) + '\n')
        ###################################################################################################
        def extract_prefix(pdb_file):
            base_name = os.path.basename(pdb_file)
            return base_name.split("FMO")[0]  

        def read_pdb_coordinates(fragment_pdb_files):
            pdb_data = {}  
            for pdb_file in fragment_pdb_files:
                with open(pdb_file, 'r') as f:
                    for line in f:
                        if line.startswith("HETATM"):  
                            cols = line.split()
                            x, y, z = map(float, [cols[5], cols[6], cols[7]])  
                            pdb_data.setdefault(pdb_file, []).append((x, y, z, pdb_file))
            return pdb_data

        def find_nearest_pdb(lig_x, lig_y, lig_z, pdb_data):
            nearest_pdb = None
            min_distance = float('inf')
            for pdb_file, coordinates in pdb_data.items():
                for x, y, z, fragment in coordinates:
                    distance = math.sqrt((x - lig_x) ** 2 + (y - lig_y) ** 2 + (z - lig_z) ** 2)
                    if distance < min_distance:
                        min_distance = distance
                        nearest_pdb = pdb_file
            return nearest_pdb, min_distance

        def find_matching_pdb(lig_x, lig_y, lig_z, pdb_data, residue, tolerance=0.1):
            for pdb_file, coordinates in pdb_data.items():
                for x, y, z, fragment in coordinates:
                    if abs(x - lig_x) <= tolerance and abs(y - lig_y) <= tolerance and abs(z - lig_z) <= tolerance:
                        return pdb_file, residue  # Exact match found
            nearest_pdb, min_distance = find_nearest_pdb(lig_x, lig_y, lig_z, pdb_data)
            if nearest_pdb:
                return nearest_pdb, residue  # Assign to the nearest PDB file
            return None, None

        def compare_files_frags(ph4_file, fragment_pdb_files, binding_energy_files, output_file):
            pdb_data = read_pdb_coordinates(fragment_pdb_files)
            binding_energies = {} 
            for energy_file in binding_energy_files:
                fragment_pdb_file = energy_file.replace("_Binding_site_energies.txt", ".pdb")  
                with open(energy_file, 'r') as f:
                    lines = f.readlines()[1:]  
                    for line in lines:
                        values = line.split()
                        residue = values[0]
                        binding_energies[(fragment_pdb_file, residue)] = values[1:]   
            with open(ph4_file, 'r') as f1, open(output_file, 'w') as output:
                lines1 = f1.readlines()
                output.write("Fragment\tTotal\tEes\tEex\tEct+mix\tEdisp\tGsol\tdistance\tAngle\tLIG.x\tLIG.y\tLIG.z\tPROT.x\tPROT.y\tPROT.z\tBond_type\tprot_atom_id\tlig_atom_id\tMatched_PDB\n")
                for line1 in lines1[1:]:  
                    values1 = line1.split()
                    lig_x, lig_y, lig_z = map(float, values1[3:6])  
                    residue = values1[0]  
                    bond_type = values1[9]  
                    matched_pdb, matched_fragment = find_matching_pdb(lig_x, lig_y, lig_z, pdb_data, residue)
                    if matched_pdb is None and "pi" in bond_type.lower():
                        matched_pdb, _ = find_nearest_pdb(lig_x, lig_y, lig_z, pdb_data)
                    energy_values = binding_energies.get((matched_pdb, residue), ["N/A"] * 6)
                    output.write("{}\t{}\t{}\t{}\n".format(
                        residue, 
                        '\t'.join(energy_values), 
                        '\t'.join(values1[1:]), 
                        matched_pdb if matched_pdb else "Not_Found"
                    ))
        if self.FE_score:
            prefix = extract_prefix(self.pdb_file)
            fragment_pdb_files_all = sorted(glob.glob(f"{prefix}FMO_lig_H_frag*_H_FE.pdb"))
            fragment_pdb_files = []
            binding_energy_files = []
            for frag_pdb in fragment_pdb_files_all:
                frag_index = os.path.basename(frag_pdb).split('_H_FE.pdb')[0]
                energy_file = f"{frag_index}_H_FE_Binding_site_energies.txt"
                if os.path.exists(energy_file):
                    fragment_pdb_files.append(frag_pdb)
                    binding_energy_files.append(energy_file)
                else:
                    print(f"⚠️ Warning: Missing energy file for {frag_pdb}, skipping.")
            # Proceed if at least one match is found
            if fragment_pdb_files and binding_energy_files:
                ph4_file_path = '../FMOPhore/Ph4_3D.txt'
                output_file_path = '../FMOPhore/Ph4_3D_FMO.txt'
                compare_files_frags(ph4_file_path, fragment_pdb_files, binding_energy_files, output_file_path)
            else:
                print("❌ No valid fragment-energy pairs found.")        

        else:
            file1_path = '../FMOPhore/Ph4_3D.txt'
            file2_path = f'{self.pdb_file[:-4]}_Binding_site_energies.txt'
            output_file_path = '../FMOPhore/Ph4_3D_FMO.txt'
            compare_files(file1_path, file2_path, output_file_path)
        ###################################################################################################
        os.chdir("../FMOPhore/")
        cwd = os.getcwd()
        data_file = 'Ph4_3D_FMO.txt'
        df_map_total = pd.read_csv(data_file, delimiter='\t', skipinitialspace=True, usecols=[0, 1])
        df_map_total = df_map_total.drop_duplicates(subset=[df_map_total.columns[0]])
        df_map_total.rename(columns={df_map_total.columns[1]: f"{self.pdb_file[:-4]}"}, inplace=True)
        selected = df_map_total[df_map_total[f"{self.pdb_file[:-4]}"] >= -400]
        selected.set_index('Fragment', inplace=True)
        selected = selected.T
        sns.set_style('darkgrid')
        plt.figure(figsize=(20, 10))
        ax = sns.heatmap(selected, linewidths=0.5, cmap='RdYlGn', center=0, fmt=".3f")
        cbar = ax.collections[0].colorbar
        cbar.set_label('kcal/mol')
        plt.tick_params(labelsize=12)
        plt.yticks(rotation='horizontal')
        plt.xlabel('')
        plt.tight_layout(pad=2)
        filename = os.path.join(cwd, 'Ph4_heatmap_lig.png')
        plt.savefig(filename)
        plt.close()
        ###################################################################################################
        def convert_ph4_to_csv(input_file='Ph4_3D_FMO.txt', output_file='Ph4_report.csv'):
            df = pd.read_csv(input_file, delim_whitespace=True)
            df = df.iloc[:, :-1]
            df.to_csv(output_file, index=False, sep=',')
        convert_ph4_to_csv()

        data_file = 'Ph4_report.csv'
        with open(data_file, 'r') as file:
            lines = file.readlines()
        lines[0] = lines[0].lstrip(',').replace(',,', ',')
        lines = [line.replace(',,', ',') for line in lines]
        with open(data_file, 'w') as file:
            file.writelines(lines)
        def format_fragment_name(fragment):
            residue = fragment[:3].capitalize()  
            chain = fragment[3] 
            number = fragment[4:]  
            return f"{residue}-{chain}-{number}"

        data_file = 'Ph4_report.csv'
        df = pd.read_csv(data_file, delimiter=',', skipinitialspace=True)
        bond_type_abbreviations = {
            'Hydrogen_Bonds': 'HB',
            'Hydrophobic_Interactions': 'Hydph',
            'Water_Bridges': 'HOH',
            'pi-Stacking': 'Pi',
            'Salt_Bridges': 'SaltB',
            'Halogen_Bonds': 'HaloB'
        }
        df['Chain_ID'] = df['Fragment'].str.extract(r'(\D+)$')[0]
        df['Residue_Number'] = df['Fragment'].str.extract(r'(\d+)')[0].astype(int)
        df_sorted = df.sort_values(by=['Chain_ID', 'Residue_Number'])
        df_grouped = df_sorted.groupby('Fragment', sort=False).agg({
            'Ees': 'first',
            'Eex': 'first',
            'Ect+mix': 'first',
            'Edisp': 'first',
            'Gsol': 'first',
            'Bond_type': lambda x: ', '.join({bond_type_abbreviations.get(item, item) for item in x})  # Collect all unique bond types and abbreviate
        }).reset_index()
        df_grouped['Formatted_Fragment'] = df_grouped['Fragment'].apply(format_fragment_name)
        fig, ax = plt.subplots(figsize=(10, 6))  # Reduced figure size to make the plot smaller
        bar_width = 0.7  
        df_grouped[['Ees', 'Eex', 'Ect+mix', 'Edisp', 'Gsol']].plot(
            kind='bar', stacked=True, color=('gold', 'mediumseagreen', 'firebrick', 'steelblue', 'darksalmon'), ax=ax, width=bar_width)
        plt.xticks(ticks=range(len(df_grouped)), labels=df_grouped['Formatted_Fragment'], rotation=90, ha='right', fontsize=18, fontweight='bold')
        plt.yticks(fontsize=18, fontweight='bold')
        # ax.set_ylim(-120, 20)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True, ncol=1, fontsize=14)
        plt.xlabel("Binding site residues", fontsize=18, fontweight='bold')
        plt.ylabel("kcal/mol", fontsize=18, fontweight='bold')
        heights = df_grouped[['Ees', 'Eex', 'Ect+mix', 'Edisp', 'Gsol']].sum(axis=1)
        for idx, (types, height) in enumerate(zip(df_grouped['Bond_type'], heights)):
            y_offset = 0.03 * height if height > 0 else 0.03 * height
            ax.text(idx, height + y_offset, types, ha='center', va='bottom' if height > 0 else 'top', fontsize=12, fontweight='bold')
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        plt.tight_layout(pad=2)
        plt.savefig("Ph4_PIEDA_adjusted.png", dpi=300) 
        plt.close()
        ################################################
        input_file = "Ph4_report.csv"
        output_file = "Ph4_lig.pdb"
        with open(input_file, "r") as file:
            lines = file.readlines()
        output_lines = []
        for i, line in enumerate(lines):
            if i == 0:  
                continue
            fields = line.strip().split(",")
            residue_id = fields[0]
            coordinates = fields[9:12]
            atom_name = fields[-1]
            residue_id = fields[0][:3]
            chain_id = fields[0][3]
            residue_number = fields[0][4:]
            output_line = "{:<6s}{:5d} {:>4s} {:<3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:>2d}\n".format(
                "ATOM", i, atom_name, residue_id, chain_id, int(residue_number), float(coordinates[0]), float(coordinates[1]),
                float(coordinates[2]), 0.00, 0.00, atom_name[0], 0)
            output_lines.append(output_line)
        with open(output_file, "w") as file:
            file.writelines(output_lines)
        ###################################################################################################

if __name__ == "__main__":
    """
    FMOPhore V 0.1 - Analysis - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    parser = argparse.ArgumentParser(description="FMOPhore Fragmention Analysis")
    parser.add_argument("-pdb", "--pdb_file", type=str, required=True, help="Path to the PDB file")
    parser.add_argument('-FE', '--FE_score', action='store_true', default=None, help=argparse.SUPPRESS)
    parser.add_argument('-BE', '--Binding_Energy', action='store_true', default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()
    pdb_processor = Analysis(args.pdb_file, args.Binding_Energy, args.FE_score)
    pdb_processor.analyze()
