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
from PIL import Image
import argparse
import pandas as pd
from tabulate import tabulate
from PIL import Image, ImageDraw, ImageFont
import math
import argparse

class Analysis:
    def __init__(self, pdb_file, Binding_Energy=False):
        self.pdb_file = pdb_file
        self.Binding_Energy = Binding_Energy
    def analyze(self):
###################################################################################################
        if self.Binding_Energy:
            with open('complex_energy.txt', 'r') as complex_file, \
                 open('protein_energy.txt', 'r') as protein_file, \
                 open('ligand_energy.txt', 'r') as ligand_file:
                complex_energy_line = complex_file.readline()
                complex_energy_match = re.search(r'-?\d+\.\d+', complex_energy_line)
                complex_energy = float(complex_energy_match.group())
                protein_energy_line = protein_file.readline()
                protein_energy_match = re.search(r'-?\d+\.\d+', protein_energy_line)
                protein_energy = float(protein_energy_match.group())
                ligand_energy_line = ligand_file.readline()
                ligand_energy_match = re.search(r'-?\d+\.\d+', ligand_energy_line)
                ligand_energy = float(ligand_energy_match.group())
            binding_energy = complex_energy - protein_energy - ligand_energy
            with open('Binding_Energy.txt', 'w') as f:
                f.write("Binding Energy (ΔE) = (E_Complex) − (E_Protein) − (E_Ligand) = {:.2f} kcal/mol\n".format(binding_energy))
        ###################################################################################################
        df = pd.read_csv('totals.csv')
        df = pd.DataFrame(df)
        df = df.T.reset_index()
        df.index.name = None
        df.columns = df.iloc[0]
        df = df.iloc[1:]
        table = tabulate(df, headers=['Energy', 'kcal/mol'])
        gray_table = table.replace('\n', '\n\033[1;30m') 
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.axis('off')
        box_width = 0.4
        box_height = 0.05
        y = 0.8  
        for i, row in enumerate(table.split('\n')):
            if 'Energy' in row:
                bbox = dict(facecolor='lightgray', boxstyle='square', pad=box_width/2)
            elif '---' in row:
                continue
            else:
                bbox = dict(facecolor='white', boxstyle='square', pad=box_width/2)
            ax.text(0.1, y-(i*box_height), row, ha='left', va='top', fontsize=12, bbox=bbox,
                    family='monospace')
        fig.savefig('Interaction_Energy.png', bbox_inches='tight')
        plt.close()
        ###################################################################################################
        with open('total_average.txt', 'r') as infile, open('DA_FMO.txt', 'w') as outfile:
            sum_delta_E = 0
            count = 0
            for line in infile:
                if "dtype: float64" in line or line.strip() == "" or "Total_average" in line or "Fragment" in line:
                    continue
                columns = line.strip().split()
                try:
                    value = float(columns[1])
                    sum_delta_E += value
                    count += 1
                except ValueError:
                    print(f"Skipping line due to conversion error: {line.strip()}")
                    continue
            if count > 0:
                avg_delta_E = sum_delta_E / count
                outfile.write("DA-FMO Interaction Energy (ΔE) = ∑(ΔE_i^FMO) / N = {:.2f} kcal/mol\n".format(avg_delta_E))
            else:
                outfile.write("No valid data to calculate the average ΔE.\n")
        ###################################################################################################
        with open('totals.csv', 'r') as infile:
            reader = csv.reader(infile)
            next(reader)  
            row = next(reader)  
            delta_E = float(row[1])  
        with open('FMO.txt', 'w') as outfile:
            outfile.write("FMO Interaction Energy (∆E) = (∆E)_ij^es + (∆E)_ij^ex + (∆E)_ij^(ct+mix) + (∆E)_ij^DI + (∆E)_^Gsol = {:.2f} kcal/mol\n".format(delta_E))
        ###################################################################################################
        if self.Binding_Energy:
            BE = []
            with open('Binding_Energy.txt', 'r') as infile:
                for line in infile:
                    if line.startswith('Binding'):
                        binding_energy_str = re.findall(r'[-+]?\d*\.\d+|\d+', line)[0]
                        binding_energy = float(binding_energy_str)
            ΔG = []
            logP = []
            with open('Desc.txt', 'r') as infile:
                for line in infile:
                    if line.startswith('logD='):
                        logD = line.strip().split('=')[1].strip()
                    elif line.startswith('logP='):
                        logP = line.strip().split('=')[1].strip()
                    elif line.startswith('nRB='):
                        nRB = float(line.strip().split('=')[1].strip())
                    ΔG = delta_E - nRB
                    ΔG = binding_energy - nRB
            with open('Affinity.txt', 'w') as outfile:
                outfile.write("Affinity = aFMO(∆E) + logP + g = a({:.2f}) + b({}) + g \n".format(delta_E, logP))
                outfile.write("ref = https://doi.org/10.1021/acsomega.2c08132 = .\n")
                outfile.write("- = -------- = -\n")
                outfile.write("- = -------- = -\n")
                outfile.write("TΔS = n(Rotatable Bonds) = {} \n".format(nRB))
                outfile.write("ΔG = (ΔE) - (TΔS) = {:.2f} kcal/mol\n".format(ΔG))
                outfile.write("ref = https://doi.org/10.1186/1758-2946-3-2 = .\n")
                outfile.write("- = -------- = -\n")
        ###################################################################################################
            with open('Energies.txt', 'w') as outfile:
                for filename in ['Binding_Energy.txt', 'FMO.txt', 'DA_FMO.txt', 'Affinity.txt']:
                # for filename in ['Binding_Energy.txt', 'FMO.txt', 'DA_FMO.txt']:
                    with open(filename, 'r') as infile:
                        outfile.write(infile.read())
                    outfile.write('- = -------- = -\n')
        ###################################################################################################
        else:
            with open('Energies.txt', 'w') as outfile:
                for filename in ['FMO.txt', 'DA_FMO.txt']:
                # for filename in ['Binding_Energy.txt', 'FMO.txt', 'DA_FMO.txt']:
                    with open(filename, 'r') as infile:
                        outfile.write(infile.read())
                    outfile.write('- = -------- = -\n')
        ###################################################################################################
        with open("Energies.txt", "r") as infile:
            lines = infile.readlines()
        lines = [line.strip() for line in lines if line.strip() != "- = -------- = -"]
        data = [line.split(" = ") for line in lines]
        df = pd.DataFrame(data, columns=["Method", "Equation","kcal/mol"])
        # print(df.to_string(index=False))
        table = df.to_string(index=False)
        table = tabulate(df, headers=["Method", "Equation","kcal/mol"])
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.axis('off')
        box_width = 0.4
        box_height = 0.05
        y = 0.8  
        for i, row in enumerate(table.split('\n')):
            if 'ref' in row:
                bbox = dict(facecolor='lightgray', boxstyle='square', pad=box_width/2)
            elif '---' in row:
                continue
            else:
                bbox = dict(facecolor='white', boxstyle='square', pad=box_width/2)
            ax.text(0.1, y-(i*box_height), row, ha='left', va='top', fontsize=12, bbox=bbox,
                    family='monospace')
        fig.savefig('Energies.png', bbox_inches='tight')
        plt.close()
        ###################################################################################################
        # create Ph4 + FMO file
        ###################################################################################################
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
        file1_path = '../FMOPhore/Ph4_3D.txt'
        file2_path = 'Binding_site_energies.txt'
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
        plt.style.use('seaborn-darkgrid')
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
        plt.show()
        ###################################################################################################
        data_file = 'Ph4_3D_FMO.txt'
        df_map_total = pd.read_csv(data_file, delimiter='\t|\s+', skipinitialspace=True)
        df_map_total.to_csv('Ph4_report.csv', index=True, sep=',')
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
            'Bond_type': lambda x: ', '.join({bond_type_abbreviations.get(item, item) for item in x})  
        }).reset_index()
        df_grouped['Formatted_Fragment'] = df_grouped['Fragment'].apply(format_fragment_name)
        fig, ax = plt.subplots(figsize=(10, 6))  
        bar_width = 0.7  
        df_grouped[['Ees', 'Eex', 'Ect+mix', 'Edisp', 'Gsol']].plot(
            kind='bar', stacked=True, color=('gold', 'mediumseagreen', 'firebrick', 'steelblue', 'darksalmon'), ax=ax, width=bar_width)
        plt.xticks(ticks=range(len(df_grouped)), labels=df_grouped['Formatted_Fragment'], rotation=90, ha='right', fontsize=18, fontweight='bold')
        plt.yticks(fontsize=18, fontweight='bold')
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
        if self.Binding_Energy:
            os.remove(f'../analysis/Binding_Energy.txt')
            os.remove(f'../analysis/complex_energy.txt')
            os.remove(f'../analysis/protein_energy.txt')
            os.remove(f'../analysis/ligand_energy.txt')
            os.remove(f'../analysis/Affinity.txt')
        os.remove(f'../analysis/PIEDA_values.txt')  
        os.remove(f'../analysis/DA_FMO.txt')
        os.remove(f'../analysis/FMO.txt') 
        os.remove(f'../analysis/total_average.txt')
        os.remove(f'Ph4_3D.txt')
if __name__ == "__main__":
    """
    FMOPhore V 0.1 - Analysis - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    pdb_processor = Analysis(self.pdb_file)
    pdb_processor.analyze()