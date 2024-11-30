###################################################################################################
# Analysis of the library 
###################################################################################################
class LibAnalysis:
    def __init__(self, pdb_file, Binding_Energy=False):
        self.pdb_file = pdb_file
        self.Binding_Energy = Binding_Energy
    def Lib_analyze(self):
        import seaborn as sns
        import pandas as pd
        import os 
        from itertools import islice
        from math import factorial
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        import re
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.offsetbox import OffsetImage, AnnotationBbox
        import csv
        from numpy import nan
        import glob
        from rdkit import Chem
        from rdkit.Chem import Draw
        import argparse
        import pandas as pd
        from tabulate import tabulate
        from PIL import Image, ImageDraw, ImageFont
        import math
        import shutil
        import numpy as np
        import matplotlib.colors as mcolors
        import matplotlib.cm as cm
        from adjustText import adjust_text

###################################################################################################
        parent_dir = os.getcwd()
        all_binding_energies = []
        if self.Binding_Energy :
            for subdir in os.listdir(parent_dir):
                subdir_path = os.path.join(parent_dir, subdir)
                if not os.path.isdir(subdir_path):
                    continue
                for subdir_chain in os.listdir(subdir_path):
                    subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                    if not os.path.isdir(subdir_chain_path):
                        continue
                    analysis_dir = os.path.join(subdir_chain_path, "analysis")
                    if os.path.isdir(analysis_dir):
                        os.chdir(analysis_dir)
                        if os.path.isfile("Energies.txt"):
                            with open('Energies.txt', 'r') as energies_file:
                                binding_energy = None
                                interaction_energy = None
                                ΔG_energy = None
                                for line in energies_file:
                                    if line.startswith('Binding'):
                                        binding_energy_str = re.findall(r'[-+]?\d*\.\d+|\d+', line)[0]
                                        binding_energy = float(binding_energy_str)
                                    elif line.startswith('FMO'):
                                        interaction_energy_str = re.findall(r'[-+]?\d*\.\d+|\d+', line)[0]
                                        interaction_energy = float(interaction_energy_str)
                                    elif line.startswith('ΔG'):
                                        ΔG_energy_str = re.findall(r'[-+]?\d*\.\d+|\d+', line)[0]
                                        ΔG_energy = float(ΔG_energy_str)
                                
                                moldescs_file = None
                                for file_name in os.listdir():
                                    if file_name.endswith(".moldescs.txt"):
                                        moldescs_file = file_name
                                        break

                                if moldescs_file:
                                    with open(moldescs_file, 'r') as moldescs_file:
                                        next(moldescs_file)
                                        desc_line = next(moldescs_file).strip().split()[1:]

                                    if all((binding_energy, interaction_energy, ΔG_energy)):
                                        all_binding_energies.append((subdir, binding_energy, interaction_energy, ΔG_energy, desc_line))
                        os.chdir(parent_dir)
            sorted_binding_energies = sorted(all_binding_energies, key=lambda x: x[1], reverse=False)

            with open('library_analysis/Data.txt', 'w') as outfile:
                outfile.write("Fragment Binding_energy Interaction_energy ΔG_energy MW Fsp3 nChiral nRings nArRings nRingSys nHBA nHBD nNEG nPOS nRB TPSA logP logD ComplexScore MolFormula smiles\n")
                for folder, binding_energy, interaction_energy, ΔG_energy, desc_line in sorted_binding_energies:
                    desc_line_str = ' '.join(desc_line).strip('[]\"')
                    outfile.write(f"{folder} {binding_energy} {interaction_energy} {ΔG_energy} {desc_line_str}\n")
        ###################################################################################################

            parent_dir = os.getcwd()

            all_binding_energies = []
            for subdir in os.listdir(parent_dir):
                subdir_path = os.path.join(parent_dir, subdir)
                if not os.path.isdir(subdir_path):
                    continue
                for subdir_chain in os.listdir(subdir_path):
                    subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                    if not os.path.isdir(subdir_chain_path):
                        continue
                    analysis_dir = os.path.join(subdir_chain_path, "analysis")
                    if os.path.isdir(analysis_dir):
                        os.chdir(analysis_dir)
                        if os.path.isfile("Energies.txt"):
                            with open('Energies.txt', 'r') as infile:
                                binding_energies = []
                                for line in infile:
                                    if line.startswith('Binding'):
                                        binding_energy_str = re.findall(r'[-+]?\d*\.\d+|\d+', line)[0]
                                        binding_energy = float(binding_energy_str)
                                        binding_energies.append(binding_energy)
                                if binding_energy is not None:
                                    all_binding_energies.append((subdir, binding_energy))
                        os.chdir(parent_dir)

            sorted_binding_energies = sorted(all_binding_energies, key=lambda x: x[1], reverse=False)

            with open('library_analysis/Binding_energies.txt', 'w') as outfile:
                outfile.write(f"Fragment   Binding_energy\n")
                for folder, energy in sorted_binding_energies:
                    outfile.write(f"{folder}    {energy}\n")

        ##################################################################################################

            all_ΔG_energies = []
            for subdir in os.listdir(parent_dir):
                subdir_path = os.path.join(parent_dir, subdir)
                if not os.path.isdir(subdir_path):
                    continue
                for subdir_chain in os.listdir(subdir_path):
                    subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                    if not os.path.isdir(subdir_chain_path):
                        continue
                    analysis_dir = os.path.join(subdir_chain_path, "analysis")
                    if os.path.isdir(analysis_dir):
                        os.chdir(analysis_dir)
                        if os.path.isfile("Energies.txt"):
                            with open('Energies.txt', 'r') as infile:
                                ΔG_energies = []
                                for line in infile:
                                    if line.startswith('ΔG'):
                                        ΔG_energy_str = re.findall(r'[-+]?\d*\.\d+|\d+', line)[0]
                                        ΔG_energy = float(ΔG_energy_str)
                                        ΔG_energies.append(ΔG_energy)
                                if ΔG_energy is not None:
                                    all_ΔG_energies.append((subdir, ΔG_energy))
                        os.chdir(parent_dir)

            sorted_ΔG_energies = sorted(all_ΔG_energies, key=lambda x: x[1], reverse=False)
            with open('library_analysis/deltaG_energies.txt', 'w') as outfile:
                outfile.write(f"Fragment   ΔG_energy\n")
                for folder, energy in sorted_ΔG_energies:
                    outfile.write(f"{folder}    {energy}\n")

        ###################################################################################################
        all_interaction_energies = []
        for subdir in os.listdir(parent_dir):
            subdir_path = os.path.join(parent_dir, subdir)
            if not os.path.isdir(subdir_path):
                continue
            for subdir_chain in os.listdir(subdir_path):
                subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                if not os.path.isdir(subdir_chain_path):
                    continue
                analysis_dir = os.path.join(subdir_chain_path, "analysis")
                if os.path.isdir(analysis_dir):
                    os.chdir(analysis_dir)
                    if os.path.isfile("Energies.txt"):
                        with open('Energies.txt', 'r') as infile:
                            interaction_energies = []
                            for line in infile:
                                if line.startswith('FMO'):
                                    interaction_energy_str = re.findall(r'[-+]?\d*\.\d+|\d+', line)[0]
                                    interaction_energy = float(interaction_energy_str)
                                    interaction_energies.append(interaction_energy)
                            if interaction_energy is not None:
                                all_interaction_energies.append((subdir_chain, interaction_energy))
                    os.chdir(parent_dir)

        sorted_interaction_energies = sorted(all_interaction_energies, key=lambda x: x[1], reverse=False)

        with open('library_analysis/Interaction_energies.txt', 'w') as outfile:
            outfile.write(f"Fragment   Interaction_energy\n")
            for folder, energy in sorted_interaction_energies:
                outfile.write(f"{folder}    {energy}\n")

        ###################################################################################################
        parent_dir = os.getcwd()
        all_bonding_energies = []
        for subdir in os.listdir(parent_dir):
            subdir_path = os.path.join(parent_dir, subdir)
            if not os.path.isdir(subdir_path):
                continue
            for subdir_chain in os.listdir(subdir_path):
                subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                if not os.path.isdir(subdir_chain_path):
                    continue
                analysis_dir = os.path.join(subdir_chain_path, "FMOPhore")
                if os.path.isdir(analysis_dir):
                    os.chdir(analysis_dir)
                    if os.path.isfile("Ph4_3D_FMO.txt"):
                        with open('Ph4_3D_FMO.txt', 'r') as infile:
                            total_sum = 0
                            duplicates = set()
                            for line in infile:
                                columns = line.strip().split('\t')
                                if len(columns) > 1:
                                    key = columns[0]  
                                    if key in duplicates:
                                        continue  
                                    try:
                                        value = float(columns[1])
                                        total_sum += value
                                        duplicates.add(key)  
                                    except ValueError:
                                        pass
                        print(f"Sum of second column {subdir}:", f"{total_sum:.2f}")
                        all_bonding_energies.append((subdir_chain, round(total_sum, 2)))
                    os.chdir(parent_dir)

        sorted_bonding_energies = sorted(all_bonding_energies, key=lambda x: x[1], reverse=False)

        with open('library_analysis/Ph4_3D_FMO.txt', 'w') as outfile:
            outfile.write(f"Fragment   Bonds_energy\n")
            for folder, energy in sorted_bonding_energies:
                outfile.write(f"{folder}    {energy}\n")

        ###################################################################################################
        # DFTB_HEATMAP_PLOT
        ###################################################################################################
        def dftb_log(parent_dir):
            os.makedirs("library_analysis", exist_ok=True)
            os.chdir("library_analysis/")
            output_file = 'ligand_heatmap.csv'
            FMO_files = glob.glob("*.log")

            selection = 'LIG'
            total_list = []
            Ees_list = []
            Eex_list = []
            Ect_list = []
            Edisp_list = []
            Gsol_list = []
            df_map_total = pd.DataFrame()  # Initialize df_map_total to ensure it exists

            for file in FMO_files:
                try:
                    pattern_01 = re.compile(r'(\s{4}I\s{4}J\sDL)')
                    pattern_02 = re.compile(r'(\s{6}NFRAG=\d+)')

                    with open(output_file, 'w+') as f:
                        for i, line in enumerate(open(file)):
                            for match in re.finditer(pattern_01, line):
                                start_line = (i + 3)

                            for match in re.finditer(pattern_02, line):
                                bla = match.groups()
                                nfrag = bla[0].split("=")[1]
                                nlines = (factorial(int(nfrag)) / (factorial(2) * factorial(int(nfrag) - 2)))

                    with open(output_file, 'w+') as f:
                        for i, line in enumerate(open(file)):
                            if i in range(start_line - 1, int(start_line + nlines)):
                                f.write(line)

                    pattern_01 = re.compile(r'(\s{6}FRGNAM\(\d\)=)')
                    pattern_02 = re.compile(r'(\s{6}INDAT)')
                    pattern_03 = re.compile(r'(\s{6})')

                    combined_pat = [
                        r'(\s{6}FRGNAM\(\d\)=)',
                        r'(\s{6}INDAT)',
                        r'(\s{6})'
                    ]
                    residues = []

                    for i, line in enumerate(open(file)):
                        for match in re.finditer(pattern_01, line):
                            start_line = i
                        for match in re.finditer(pattern_02, line):
                            end_line = i

                    for i, line in enumerate(open(file)):
                        if i in range(start_line, end_line):
                            residues.append(re.sub('|'.join(combined_pat), " ", line))

                    residues = list(map(lambda x: x.strip(), residues))
                    residues = [x.strip(' ') for x in residues]

                    my_residues = []
                    for element in residues:
                        my_residues.extend(element.split(','))

                    my_residues = [x.strip(' ') for x in my_residues]
                    my_residues = [x for x in my_residues if x != '']

                    fragment_num = str(my_residues.index(selection) + 1)

                    df = pd.DataFrame()
                    for i, line in enumerate(open(output_file)):
                        if fragment_num in (line[:10]):
                            temp = pd.DataFrame({
                                "I": [line[:5]],
                                "J": [line[6:10]],
                                "DL": [line[11:13]],
                                "Z": [line[14:16]],
                                "R": [line[17:23]],
                                "Q(I->J)": [line[24:31]],
                                "EIJ-EI-EJ": [line[32:41]],
                                "dDIJ*VIJ": [line[42:50]],
                                "total": [line[51:60]],
                                "Ees": [line[61:70]],
                                "Eex": [line[71:79]],
                                "Ect+mix": [line[80:88]],
                                "Edisp": [line[89:97]],
                                "Gsol": [line[98:106]]
                            })
                            df = pd.concat([df, temp])

                    my_residues.pop(int(fragment_num) - 1)
                    df.insert(0, "Fragment", my_residues)

                    df_map = df[['Fragment', 'total']].T
                    df_map.columns = df_map.iloc[0]
                    df_map = df_map.drop(['Fragment'], axis=0)
                    df_map = df_map.set_index([pd.Index([os.path.splitext(file)[0]])])
                    df_map = df_map.apply(pd.to_numeric, errors='coerce')

                    if df_map_total.empty:
                        df_map_total = df_map
                    else:
                        df_map_total = df_map_total.append(df_map, sort=False)
                        df_map_total = df_map_total.reindex(sorted(df_map_total.columns, key=lambda s: s[s.rfind('-'):]), axis=1)
                        df_map_total = df_map_total.reindex(sorted(df_map_total.columns, key=lambda x: int(re.search(r'\d+$', x).group())), axis=1)

                        df_map_total.index = df_map_total.index.str.replace('_FMO', '')
                        df_map_total.index = df_map_total.index.str.replace('min', '')

                    cols = df.columns.drop(['Fragment', 'DL'])
                    df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')

                    df2 = df

                    sorter = list(df_map_total.columns)
                    true_sort = [s for s in sorter if s in df2.Fragment.unique()]
                    df2 = df2.set_index('Fragment').loc[true_sort].reset_index()
                    df2.set_index("Fragment", drop=True, inplace=True)

                    if 'Ees' not in locals():
                        Ees = [list(df2['Ees'])]
                        Eex = [list(df2['Eex'])]
                        Ect = [list(df2['Ect+mix'])]
                        Edisp = [list(df2['Edisp'])]
                        Gsol = [list(df2['Gsol'])]
                    else:
                        Ees.append(list(df2['Ees']))
                        Eex.append(list(df2['Eex']))
                        Ect.append(list(df2['Ect+mix']))
                        Edisp.append(list(df2['Edisp']))
                        Gsol.append(list(df2['Gsol']))

                    total_list.append(df['total'].sum())
                    Ees_list.append(df['Ees'].sum())
                    Eex_list.append(df['Eex'].sum())
                    Ect_list.append(df['Ect+mix'].sum())
                    Edisp_list.append(df['Edisp'].sum())
                    Gsol_list.append(df['Gsol'].sum())

                except Exception as e:
                    print(f"Error processing file {file}: {e}")
                    continue

            if not df_map_total.empty:
                bla = list(map(tuple, np.where(np.isnan(df_map_total.values))))
                bla_a = bla[0]
                bla_b = bla[1]

                for f, b in zip(bla_a, bla_b):
                    Ees[f].insert(b, nan)
                    Eex[f].insert(b, nan)
                    Ect[f].insert(b, nan)
                    Edisp[f].insert(b, nan)

                df_map_total.insert(0, 'total', total_list)
                df_map_total = df_map_total.sort_values(by=['total'])

                df_map_total.to_csv('export_dataframe.csv', header=True)

                selected = df_map_total[df_map_total["total"] >= -400]
                selected = selected.drop(columns=['total'])

                plt.style.use('seaborn-darkgrid')
                plt.figure(figsize=(20, 10))
                ax = sns.heatmap(selected, linewidths=.5, cmap="RdYlGn", center=0, fmt=".0f")
                cbar = ax.collections[0].colorbar
                cbar.set_label("kcal/mol")

                plt.tick_params(labelsize=12)
                plt.yticks(rotation='horizontal')
                plt.xlabel('')
                plt.tight_layout(pad=2)
                plt.savefig('DFTB_heatmap_lig.png')

                print("Heatmap has been saved as 'DFTB_heatmap_lig.png'")

        parent_dir = "./"
        dftb_log(parent_dir)
        ##############################################################################
        os.chdir("..")
        # Function to format fragment names
        def format_fragment_name(fragment):
            # Split the fragment into parts: residue (first 3 letters), chain (letter), and number (numeric part)
            residue = fragment[:3].capitalize()
            chain = fragment[3]
            number = fragment[4:]
            # Return the formatted string in the form Residue-Chain-Number
            return f"{residue}-{chain}-{number}"

        parent_dir = os.getcwd()
        all_selected_data = pd.DataFrame()

        # Loop through the directories and process data
        for subdir in os.listdir(parent_dir):
            subdir_path = os.path.join(parent_dir, subdir)
            if not os.path.isdir(subdir_path):
                continue
            for subdir_chain in os.listdir(subdir_path):
                subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                if not os.path.isdir(subdir_chain_path):
                    continue
                analysis_dir = os.path.join(subdir_chain_path, "FMOPhore")
                if os.path.isdir(analysis_dir):
                    os.chdir(analysis_dir)
                    for file in os.listdir(analysis_dir):
                        if file == 'Ph4_report.csv':
                            file_path = os.path.join(analysis_dir, file)
                            df_map_total = pd.read_csv(file_path, delimiter=',', skipinitialspace=True, usecols=[0, 1])
                            df_map_total = df_map_total.drop_duplicates(subset=[df_map_total.columns[0]])
                            df_map_total.rename(columns={df_map_total.columns[1]: subdir_chain}, inplace=True)
                            selected = df_map_total[df_map_total[subdir_chain] >= -400]
                            selected.set_index('Fragment', inplace=True)
                            selected = selected.T
                            all_selected_data = all_selected_data.append(selected)
                os.chdir(parent_dir)

        library_analysis_path = os.path.join(parent_dir, 'library_analysis')
        if not os.path.exists(library_analysis_path):
            os.makedirs(library_analysis_path)
            
        output_filename = os.path.join(parent_dir, library_analysis_path, 'selected_data.csv')
        all_selected_data.to_csv(output_filename)
        selected_data_path = os.path.join(parent_dir, library_analysis_path, 'selected_data.csv')
        selected_data = pd.read_csv(selected_data_path, index_col=0)
        selected_data['total'] = selected_data.sum(axis=1)
        selected_data = selected_data.sort_values('total', ascending=True)
        selected_data = selected_data.drop('total', axis=1)

        def get_chain_id_and_residue_number(column_name):
            chain_id = column_name[3]
            residue_number = column_name[4:]
            return (chain_id, int(residue_number))

        sorted_columns = sorted(selected_data.columns, key=get_chain_id_and_residue_number)
        selected_data = selected_data[sorted_columns]
        percentage_found = (selected_data.count() / len(selected_data) * 100).round(0).astype(int)

        # Apply the format_fragment_name function to the sorted_columns to format x-axis labels
        formatted_columns = [format_fragment_name(col) for col in sorted_columns]

        plt.style.use('seaborn-darkgrid')
        plt.figure(figsize=(20, 10))

        # Create the heatmap
        ax = sns.heatmap(selected_data,
                         linewidths=0.5,
                         linecolor='grey',
                         cmap='coolwarm',
                         center=0, 
                         annot=True,
                         annot_kws={"size": 8},
                         fmt=".1f",
                         vmax=30, 
                         vmin=-50)

        # Draw vertical lines to separate the columns
        ax.vlines(range(len(selected_data.columns) + 1), *ax.get_ylim(), colors='grey', linewidth=0.5)
        ax.set_facecolor('white')

        # Add a second x-axis to display percentage found for each binding site
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        tick_positions = np.arange(0.5, len(selected_data.columns), 1)
        ax2.set_xticks(tick_positions)
        ax2.set_xticks(range(len(percentage_found)))
        ax2.set_xticklabels(percentage_found, ha='center', rotation=0, fontsize=20) 

        # Set colorbar label and other aesthetics
        cbar = ax.collections[0].colorbar
        cbar.set_label('kcal/mol')

        # Increase tick label size for both axes
        plt.tick_params(labelsize=15)

        # Rotate y-axis labels and set font size
        plt.yticks(rotation=0, fontsize=20)  # Set y-axis labels to horizontal

        # Increase font size for the x and y axis labels
        plt.xlabel('Binding Site Residues', fontsize=20)  # Increase x-axis label font size
        plt.ylabel('PDB IDs', fontsize=20)  # Increase y-axis label font size

        # Set the title and labels for the heatmap with increased font size
        plt.title('Percentage of Interaction %', fontsize=20)  # Increase title font size
        plt.xticks(rotation=0, ha="center", fontsize=20)  # Increase x-axis tick label size
        plt.yticks(rotation=0, fontsize=20)  # Increase y-axis tick label size

        # Format x-axis labels and rotate them to 90 degrees
        ax.set_xticklabels(formatted_columns, fontweight='bold', rotation=90, fontsize=20)

        # Set bold font for the y-axis labels
        ax.set_yticklabels(ax.get_yticklabels(), fontweight='bold', fontsize=20)

        # Adjust layout to ensure no labels are cut off
        plt.tight_layout(pad=3)  # Adjust padding for better layout

        # Save the heatmap
        heatmap_filename = os.path.join(parent_dir, library_analysis_path, 'Ph4_heatmap.png')
        plt.savefig(heatmap_filename, bbox_inches='tight')  # Save with tight layout to avoid cutting off

        plt.savefig(heatmap_filename)
        ###############################################################################
        library_analysis_path = os.path.join(parent_dir, 'library_analysis')
        output_filename = os.path.join(library_analysis_path, 'selected_data.csv')
        df = pd.read_csv(output_filename, index_col=0)
        def get_chain_id_and_residue_number(column_name):
            chain_id = column_name[3]
            residue_number = column_name[4:]
            return (chain_id, int(residue_number))

        sorted_columns = sorted(df.columns, key=get_chain_id_and_residue_number)
        df = df[sorted_columns]
        percentage_values = df.notna().mean() * 100
        lowest_energy = df.min()
        sorted_percentage_values = percentage_values.reindex(sorted_columns)
        sorted_lowest_energy = lowest_energy.reindex(sorted_columns)
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

        bars = ax.bar(sorted_percentage_values.index, sorted_percentage_values, color=plt.cm.coolwarm((sorted_lowest_energy - (-50)) / (30 - (-50))), edgecolor='black')
        for bar, percentage, energy in zip(bars, sorted_percentage_values, sorted_lowest_energy):
            height = bar.get_height()
            ax.annotate(f'{percentage:.1f}%', xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords="offset points", ha='center', va='bottom')
            ax.annotate(f'{energy:.2f}', xy=(bar.get_x() + bar.get_width() / 2, height / 2),
                        xytext=(0, 0), textcoords="offset points", ha='center', va='center', color='white')

        sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=-50, vmax=30))
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.set_label('Lowest Energy Total')
        ax.set_xlabel('Binding Site Residues')
        ax.set_ylabel('Percentage of interaction (%)')
        ax.set_title('Binding site hotspots')
        plt.xticks(rotation=90)
        plt.tight_layout()
        heatmap_filename = os.path.join(library_analysis_path, 'Binding_site_residues_plot.png')
        plt.savefig(heatmap_filename, facecolor=fig.get_facecolor())

        ###############################################################################
        def keep_lowest_rows(file_path):
            with open(file_path, 'r') as file:
                reader = csv.reader(file, delimiter=',')
                data = list(reader)

            header = data[0]

            lowest_values = {}
            for i, column in enumerate(header[1:]):
                column_values = [float(row[i+1]) for row in data[1:] if row[i+1]]
                if column_values:
                    lowest_value = min(column_values)
                    lowest_values[column] = lowest_value

            new_data = [header]
            for row in data[1:]:
                new_row = [row[0]]
                for i, value in enumerate(row[1:]):
                    column = header[i+1]
                    if column in lowest_values and value and float(value) == lowest_values[column]:
                        new_row.append(value)
                    else:
                        new_row.append('')
                new_data.append(new_row)

            output_path = os.path.join(parent_dir, library_analysis_path, 'lowest_bondE.csv')
            with open(output_path, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(new_data)

        output_filename = os.path.join(parent_dir, library_analysis_path, 'selected_data.csv')
        keep_lowest_rows(output_filename)
        ###############################################################################
        library_analysis_path = os.path.join(parent_dir, 'library_analysis')
        lowest_row_data_path = os.path.join(parent_dir, library_analysis_path, 'lowest_bondE.csv')
        selected_data = pd.read_csv(lowest_row_data_path, index_col=0)
        selected_data['total'] = selected_data.sum(axis=1)
        selected_data = selected_data[selected_data['total'] < 0] 
        selected_data = selected_data.sort_values('total', ascending=True)
        selected_data = selected_data.drop('total', axis=1)

        def get_chain_id_and_residue_number(column_name):
            chain_id = column_name[3]
            residue_number = column_name[4:]
            return (chain_id, int(residue_number))

        sorted_columns = sorted(selected_data.columns, key=get_chain_id_and_residue_number)
        selected_data = selected_data[sorted_columns]

        plt.style.use('seaborn-darkgrid')
        plt.figure(figsize=(20, 10))
        ax = sns.heatmap(selected_data, linewidths=0.5, cmap='RdYlGn', center=0, fmt=".3f")
        cbar = ax.collections[0].colorbar
        cbar.set_label('kcal/mol')
        plt.tick_params(labelsize=12)
        plt.yticks(rotation='horizontal')
        plt.xlabel('')
        plt.tight_layout(pad=2)

        heatmap_filename = os.path.join(parent_dir, library_analysis_path, 'Ph4_heatmap_suggested.png')
        plt.savefig(heatmap_filename)

        output_folder = os.path.join(parent_dir, library_analysis_path, 'suggested_ligands')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        with open(lowest_row_data_path, 'r') as file:
            next(file)
            for line in file:
                fields = line.strip().split(',')
                if any(field.strip() and float(field) < 0 for field in fields[1:]):
                    folder_prefix = fields[0]
                    print(folder_prefix)
                    for root, dirs, files in os.walk('.'):
                        for dir_name in dirs:
                            if dir_name.startswith(folder_prefix):
                                full_dir_path = os.path.join(root, dir_name)
                                print(full_dir_path)
                                for filename in os.listdir(full_dir_path):
                                    if filename.startswith(folder_prefix) and filename.endswith('_lig_H.pdb'):
                                        shutil.copy(os.path.join(full_dir_path, filename), output_folder)
                                    if filename.endswith('Ph4_PIEDA.png'):
                                        new_filename = f"{folder_prefix}_Ph4_PIEDA.png"
                                        shutil.copy(os.path.join(full_dir_path, filename), os.path.join(output_folder, new_filename))
                                    if filename.endswith('Ph4_3D_FMO.txt'):
                                        new_filename = f"{folder_prefix}_Ph4_3D_FMO.txt"
                                        shutil.copy(os.path.join(full_dir_path, filename), os.path.join(output_folder, new_filename))
        ###############################################################################
        # merge Ph4s
        # ###############################################################################
        # parent_dir = os.getcwd()
        # output_file = "library_analysis/Ph4_reports.csv"
        # all_data = []
        # error_file = "errors.txt"
        # all_selected_data = pd.DataFrame()
        # for subdir in os.listdir(parent_dir):
        #     subdir_path = os.path.join(parent_dir, subdir)
        #     if not os.path.isdir(subdir_path):
        #         continue
        #     for subdir_chain in os.listdir(subdir_path):
        #         subdir_chain_path = os.path.join(subdir_path, subdir_chain)
        #         if not os.path.isdir(subdir_chain_path):
        #             continue
        #         file_path = os.path.join(subdir_chain_path, "FMOPhore", "Ph4_report.csv")
        #         if os.path.isfile(file_path):
        #             with open(file_path, "r") as file:
        #                 csv_reader = csv.reader(file)
        #                 next(csv_reader)  
        #                 folder_name = subdir  
        #                 for row in csv_reader:
        #                     all_data.append(row + [folder_name])

        def process_ph4_reports(parent_dir, error_file, all_data):
            for subdir in os.listdir(parent_dir):
                subdir_path = os.path.join(parent_dir, subdir)
                if not os.path.isdir(subdir_path):
                    continue
                for subdir_chain in os.listdir(subdir_path):
                    subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                    if not os.path.isdir(subdir_chain_path):
                        continue
                    file_path = os.path.join(subdir_chain_path, "FMOPhore", "Ph4_report.csv")
                    if not os.path.isfile(file_path):
                        with open(error_file, "a") as error_f:
                            error_f.write(f"{subdir_chain_path}\n")
                        continue
                    with open(file_path, "r") as file:
                        csv_reader = csv.reader(file)
                        next(csv_reader) 
                        folder_name = subdir
                        for row in csv_reader:
                            all_data.append(row + [folder_name])

        parent_dir = os.getcwd() 
        print(parent_dir)
        error_file = "errors_library.log"
        all_data = []
        process_ph4_reports(parent_dir, error_file, all_data)
        output_file = "library_analysis/Ph4_reports.csv"
        # os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, "w", newline="") as file:
            csv_writer = csv.writer(file)
            csv_writer.writerow(["Fragment", "total", "Ees", "Eex", "Ect+mix", "Edisp", "Gsol",
                                 "distance", "Angle", "LIG.x", "LIG.y", "LIG.z", "PROT.x",
                                 "PROT.y", "PROT.z", "Bond_type", "prot_atom_id", "lig_atom_id", "PDB_complex"])
            csv_writer.writerows(all_data)
        print("Data from all files has been combined and saved to Ph4_reports.csv.")
        ###############################################################################
        # Filter Ph4s
        ###############################################################################
        input_file = "library_analysis/Ph4_reports.csv"
        output_file = "library_analysis/Ph4_reports_filtered.csv"

        lowest_totals = {}

        with open(input_file, "r") as file:
            csv_reader = csv.reader(file)
            header = next(csv_reader)  

            for row in csv_reader:
                residue_id = row[0]  
                bond_type = row[-4]  
                total = float(row[1])  

                if (residue_id, bond_type) in lowest_totals:
                    if total < lowest_totals[(residue_id, bond_type)]:
                        lowest_totals[(residue_id, bond_type)] = total
                else:
                    lowest_totals[(residue_id, bond_type)] = total

        filtered_rows = []
        with open(input_file, "r") as file:
            csv_reader = csv.reader(file)
            next(csv_reader)  
            for row in csv_reader:
                residue_id = row[0]  
                bond_type = row[-4]  
                total = float(row[1])  

                if total == lowest_totals[(residue_id, bond_type)]:
                    filtered_rows.append(row)

        # sorted_rows = sorted(filtered_rows, key=lambda x: int(x[0].split("-")[1]))
        sorted_rows = filtered_rows

        with open(output_file, "w", newline="") as file:
            csv_writer = csv.writer(file)
            csv_writer.writerow(header)
            csv_writer.writerows(sorted_rows)

        print("Filtered data has been saved to filtered_Ph4_reports.csv.")
        ##############################################################################
        input_file = "library_analysis/Ph4_reports_filtered.csv"
        output_file = "library_analysis/Ph4_lig.pdb"

        with open(input_file, "r") as file:
            lines = file.readlines()

        output_lines = []
        for i, line in enumerate(lines):
            if i == 0:  
                continue
            
            fields = line.strip().split(",")
            residue_id = fields[0]
            coordinates = fields[9:12]
            atom_name = fields[-2]
            residue_number = fields[0][4:]
            residue_id = fields[0][:3]
            chain_id = fields[0][3:4]
            # residue_id = "LIG"
            # print(chain_id)
            output_line = "{:<6s}{:5d} {:>4s} {:<3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:>2d}\n".format(
            "ATOM", i, atom_name, residue_id, chain_id, int(residue_number), float(coordinates[0]), float(coordinates[1]),float(coordinates[2]), 0.00, 0.00, atom_name[0], 0)
            output_lines.append(output_line)

        with open(output_file, "w") as file:
            file.writelines(output_lines)

        print("File converted successfully.")
        ###############################################################################
        input_file = "library_analysis/Ph4_reports.csv"
        output_file = "library_analysis/Ph4_lig_cloud.pdb"

        with open(input_file, "r") as file:
            lines = file.readlines()

        output_lines = []

        for i, line in enumerate(lines):
            if i == 0:  
                continue
            
            fields = line.strip().split(",")
            residue_id = fields[0]
            coordinates = fields[9:12]
            atom_name = fields[-2]
            residue_number = fields[0][4:]
            residue_id = fields[0][:3]
            chain_id = fields[0][3:4]

            output_line = "{:<6s}{:5d} {:>4s} {:<3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:>2d}\n".format(
                "ATOM", i, atom_name, residue_id, chain_id, int(residue_number), float(coordinates[0]), float(coordinates[1]),
                float(coordinates[2]), 0.00, 0.00, atom_name[0], 0)
            output_lines.append(output_line)

        with open(output_file, "w") as file:
            file.writelines(output_lines)
        ###############################################################################
        input_file = "library_analysis/Ph4_reports.csv"
        output_file_BS = "library_analysis/Ph4_BS_Dynamics.pdb"

        with open(input_file, "r") as file:
            lines = file.readlines()

        output_lines_BS = []

        for i, line in enumerate(lines):
            if i == 0:  
                continue

            fields = line.strip().split(",")  
            residue_id_BS = fields[0]
            coordinates_BS = fields[12:15]
            atom_name_BS = fields[-3]
            chain_id_BS = fields[0][3:4]
            residue_number_BS = fields[0][4:]
            residue_id_BS = fields[0][:3]

            output_line_BS = "{:<6s}{:5d} {:>4s} {:<3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:>2d}\n".format(
                "ATOM", i, atom_name_BS, residue_id_BS, chain_id_BS, int(residue_number_BS), float(coordinates_BS[0]), 
                float(coordinates_BS[1]), float(coordinates_BS[2]), 0.00, 0.00, atom_name_BS[0], 0)
            output_lines_BS.append(output_line_BS)

        with open(output_file_BS, "w") as file:
            file.writelines(output_lines_BS)

        print("File converted successfully.")

        ##############################################################################
        first_file_path = 'library_analysis/Ph4_reports_filtered.csv'
        residue_ids_numbers = []
        with open(first_file_path, 'r') as file:
            lines = file.readlines()
            for line in lines[1:]: 
                residue_id_number = line.split(',')[0].strip()
                residue_id, chain_id, residue_number = residue_id_number[:3], residue_id_number[3], residue_id_number[4:]
                residue_ids_numbers.append((residue_id, chain_id, residue_number))

        directory_path = './'
        extracted_residues = []
        for root, dirs, files in os.walk(directory_path):
            for file in files:
                if file.endswith('_H.pdb') and not file.endswith('lig_H.pdb'):
                    file_path = os.path.join(root, file)
                    if not os.path.isfile(file_path):  # Check if the file exists
                        continue
                    with open(file_path, 'r') as pdb_file:
                        pdb_lines = pdb_file.readlines()
                        for line in pdb_lines:
                            if line.startswith('ATOM'):
                                residue_id = line[17:20].strip()
                                chain_id = line[21:22].strip()
                                residue_number = line[22:26].strip()
                                if (residue_id, chain_id, residue_number) in residue_ids_numbers:
                                    extracted_residues.append(line.strip())

        output_file_path = 'library_analysis/BS_dynamics.pdb'
        with open(output_file_path, 'w') as output_file:
            for residue_line in extracted_residues:
                output_file.write(residue_line + '\n')

        ###############################################################################
        output_file = "library_analysis/pocket.pdb"
        written_lines = set() 
        with open(output_file, "w") as output_f:
            for subdir in os.listdir(parent_dir):
                subdir_path = os.path.join(parent_dir, subdir)
                if not os.path.isdir(subdir_path):
                    continue
                for subdir_chain in os.listdir(subdir_path):
                    subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                    if not os.path.isdir(subdir_chain_path):
                        continue
                    analysis_BS_path = os.path.join(subdir_chain_path, "FMOPhore")
                    if not os.path.isdir(analysis_BS_path):
                        continue 
                    if os.path.isdir(analysis_BS_path):
                        for filename in os.listdir(analysis_BS_path):
                            if filename == "pocket.pdb":
                                file_path = os.path.join(analysis_BS_path, filename)
                                with open(file_path, "r") as input_f:
                                    lines1 = input_f.readlines()
                                    lines1 = [line for line in lines1 if 'HEADER' not in line]
                                    lines1 = [line for line in lines1 if 'END' not in line]
                                    lines1 = [line for line in lines1 if 'TER' not in line]
                                    for line in lines1:
                                        if line not in written_lines:
                                            output_f.write(line)
                                            written_lines.add(line)
        output_file = "library_analysis/pocket.pqr"
        written_lines = set() 
        with open(output_file, "w") as output_f:
            for subdir in os.listdir(parent_dir):
                subdir_path = os.path.join(parent_dir, subdir)
                if not os.path.isdir(subdir_path):
                    continue
                for subdir_chain in os.listdir(subdir_path):
                    subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                    if not os.path.isdir(subdir_chain_path):
                        continue
                    analysis_BS_path = os.path.join(subdir_chain_path, "FMOPhore")
                    if not os.path.isdir(analysis_BS_path):
                        continue 
                    if os.path.isdir(analysis_BS_path):
                        for filename in os.listdir(analysis_BS_path):
                            if filename == "pocket.pqr":
                                file_path = os.path.join(analysis_BS_path, filename)
                                with open(file_path, "r") as input_f:
                                    lines1 = input_f.readlines()
                                    lines1 = [line for line in lines1 if 'HEADER' not in line]
                                    lines1 = [line for line in lines1 if 'END' not in line]
                                    lines1 = [line for line in lines1 if 'TER' not in line]
                                    for line in lines1:
                                        if line not in written_lines:
                                            output_f.write(line)
                                            written_lines.add(line)

        ###############################################################################
        def distance(coord1, coord2):
            return np.linalg.norm(np.array(coord2) - np.array(coord1))

        def process_directory_structure(parent_dir):
            for subdir in os.listdir(parent_dir):
                subdir_path = os.path.join(parent_dir, subdir)
                if not os.path.isdir(subdir_path):
                    return
                process_subdirectories(subdir_path)

        def process_subdirectories(subdir_path):
            for subdir_chain in os.listdir(subdir_path):
                subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                if not os.path.isdir(subdir_chain_path):
                    return
                process_files(subdir_chain_path)

        def process_files(subdir_chain_path):
            analysis_dir = os.path.join(subdir_chain_path, "analysis")
            ligand_files = glob.glob(os.path.join(subdir_chain_path, "*_FMO_lig_H.pdb"))
            if not os.path.isdir(analysis_dir):
                return
            for ligand_file in ligand_files:
                if "FMOPhore" not in ligand_file:
                    matching_ligand = ligand_file
                    process_ligand_file(matching_ligand, analysis_dir)
                    break

        def process_ligand_file(matching_ligand, analysis_dir):
            with open(matching_ligand, 'r') as file:
                file2_coords = np.array([np.fromstring(line[30:54], sep=' ') for line in file if line.startswith('HETATM')])

            all_growth_vectors_path = os.path.join(analysis_dir, '../../../library_analysis/pocket.pqr')
            pocket_path = os.path.join(analysis_dir, '../../../library_analysis/pocket.pdb')
            single_growth_vectors_path = os.path.join(analysis_dir, 'growth_vectors.pdb')

            with open(pocket_path, 'r') as file:
                pocket_coords = np.array([np.fromstring(line[30:54], sep=' ') for line in file if line.startswith('ATOM')])

            process_vectors(file2_coords, all_growth_vectors_path, single_growth_vectors_path)

        def process_vectors(file2_coords, all_growth_vectors_path, single_growth_vectors_path):
            with open(all_growth_vectors_path, 'r') as file:
                growth_vectors = [line for line in file if line.startswith('ATOM')]
            
            filtered_coords = []
            for line in growth_vectors:
                coords = np.fromstring(line[30:54], sep=' ')
                if np.all(np.linalg.norm(file2_coords - coords, axis=1) > 2):  # Using broadcasting for distance check
                    filtered_coords.append(line)

            with open(single_growth_vectors_path, 'w') as output_file:
                output_file.writelines(filtered_coords)

        parent_dir = os.getcwd()
        process_directory_structure(parent_dir)
        ###############################################################################

        ##########################################################################################################################
        color_mapping = {
            'pi-Stacking': 'orange',
            'Hydrophobic_Interactions': 'green',
            'Hydrogen_Bonds': 'blue',
            'Halogen_Bonds': 'cyan',
            'Salt_Bridges': 'red',
            'Water_Bridges': 'magenta',
            'pi-Cation_Interactions': 'yellow',
            'Metal_Complexes': 'brown'
        }

        file_path = "library_analysis/Ph4_reports.csv"
        df = pd.read_csv(file_path)
        def get_chain_id_and_residue_number(fragment_name):
            chain_id = fragment_name[3]
            residue_number = fragment_name[4:]
            return (chain_id, int(residue_number))
        sorted_fragments = sorted(df['Fragment'].unique(), key=get_chain_id_and_residue_number)
        fragment_occurrence_counts = df['Fragment'].value_counts()
        fragment_occurrence_percentages = (fragment_occurrence_counts / len(df)) * 100
        bond_type_counts = df.groupby(['Fragment', 'Bond_type']).size().unstack(fill_value=0)
        bond_type_percentages = bond_type_counts.div(bond_type_counts.sum(axis=1), axis=0)
        fragment_occurrence_percentages = fragment_occurrence_percentages.reindex(sorted_fragments)
        bond_type_percentages = bond_type_percentages.reindex(sorted_fragments)
        fig, ax = plt.subplots(figsize=(14, 8))
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

        bottom = np.zeros(len(sorted_fragments))
        for bond_type in bond_type_percentages.columns:
            color = color_mapping.get(bond_type, 'gray')
            heights = fragment_occurrence_percentages * bond_type_percentages[bond_type]
            ax.bar(sorted_fragments, heights, label=bond_type, color=color, edgecolor='black', bottom=bottom)
            bottom += heights

        for i, fragment in enumerate(sorted_fragments):
            total_percentage = fragment_occurrence_percentages[fragment]
            ax.annotate(f'{total_percentage:.1f}%', xy=(i, total_percentage),
                        xytext=(0, 3), textcoords="offset points", ha='center', va='bottom')

        ax.legend(title='Bond Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12, title_fontsize=14)
        ax.set_xlabel('Fragment', fontsize=14, fontweight='bold')
        ax.set_ylabel('Percentage of Occurrence (%)', fontsize=14, fontweight='bold')
        ax.set_title('Percentage of Occurrence and Bond Types for Each Fragment', fontsize=16, fontweight='bold')
        plt.xticks(rotation=90, fontsize=12, fontweight='bold')
        plt.yticks(fontsize=12, fontweight='bold')
        plt.tight_layout()
        output_file = 'library_analysis/Binding_site_residues_plot_with_bond_types_and_occurrence.png'
        plt.savefig(output_file)
        output_data = bond_type_percentages.multiply(fragment_occurrence_percentages, axis=0)
        output_data.insert(0, 'Fragment', fragment_occurrence_percentages.index)
        output_data.reset_index(drop=True, inplace=True)
        output_csv_file = 'library_analysis/binding_site_residues_bond_type_percentages.csv'
        output_data.to_csv(output_csv_file, index=False, float_format='%.1f')
        ##########################################################################################################################
        library_analysis_path = os.path.join(parent_dir, 'library_analysis')
        file_path = os.path.join(library_analysis_path, 'selected_data.csv')
        df = pd.read_csv(file_path, index_col=0)
        def get_chain_id_and_residue_number(column_name):
            chain_id = column_name[3]
            residue_number = column_name[4:]
            return (chain_id, int(residue_number))

        sorted_columns = sorted(df.columns, key=get_chain_id_and_residue_number)
        df = df[sorted_columns]
        percentage_values = df.notna().mean() * 100
        lowest_energy = df.min()
        sorted_percentage_values = percentage_values.reindex(sorted_columns)
        sorted_lowest_energy = lowest_energy.reindex(sorted_columns)
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

        bars = ax.bar(sorted_percentage_values.index, sorted_percentage_values, color=plt.cm.coolwarm((sorted_lowest_energy - (-50)) / (30 - (-50))), edgecolor='black')
        for bar, percentage, energy in zip(bars, sorted_percentage_values, sorted_lowest_energy):
            height = bar.get_height()
            ax.annotate(f'{percentage:.1f}%', xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=12, fontweight='bold')
            ax.annotate(f'{energy:.2f}', xy=(bar.get_x() + bar.get_width() / 2, height / 2),
                        xytext=(0, 0), textcoords="offset points", ha='center', va='center', color='white', fontsize=12, fontweight='bold')

        sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=-50, vmax=30))
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.set_label('Lowest Energy Total', fontsize=14, fontweight='bold')
        ax.set_xlabel('Binding Site Residues', fontsize=14, fontweight='bold')
        ax.set_ylabel('Percentage of interaction (%)', fontsize=14, fontweight='bold')
        ax.set_title('Binding site hotspots', fontsize=16, fontweight='bold')
        plt.xticks(rotation=90, fontsize=12, fontweight='bold')
        plt.yticks(fontsize=12, fontweight='bold')  # Ensure y-ticks are bold and larger
        plt.tight_layout()
        # output_file_coolwarm = 'library_analysis/Binding_site_residues_plot1.png'
        # plt.savefig(output_file_coolwarm)
        output_data = pd.DataFrame({
            'Binding Site Residue': sorted_percentage_values.index,
            'Percentage of Interaction (%)': sorted_percentage_values.values,
            'Lowest Energy': sorted_lowest_energy.values
        })
        output_csv_file = 'library_analysis/binding_site_residues_interaction_and_energy.csv'
        output_data.to_csv(output_csv_file, index=False, float_format='%.1f')
        # print(f"Plot saved as {output_file_coolwarm}")
        print(f"Data saved as {output_csv_file}")
        ##########################################################################################################################
        file1_path = 'library_analysis/binding_site_residues_bond_type_percentages.csv'
        df1 = pd.read_csv(file1_path)
        file2_path = 'library_analysis/binding_site_residues_interaction_and_energy.csv'
        df2 = pd.read_csv(file2_path)
        merged_df = pd.merge(df1, df2, how='inner', left_on='Fragment', right_on='Binding Site Residue')
        merged_df.drop(columns=['Binding Site Residue'], inplace=True)
        sorted_fragments = sorted(merged_df['Fragment'].unique(), key=get_chain_id_and_residue_number)
        merged_df = merged_df.set_index('Fragment').reindex(sorted_fragments).reset_index()
        output_file = 'library_analysis/merged_file.csv'
        merged_df.to_csv(output_file, index=False)

        print(f"Merged data saved as {output_file}")


        def format_fragment_name(fragment):
            residue = fragment[:3].capitalize()
            chain = fragment[3]
            number = fragment[4:]
            return f"{residue}-{chain}-{number}"

        data = pd.read_csv('library_analysis/merged_file.csv')
        min_percentage = data['Percentage of Interaction (%)'].min()
        max_percentage = data['Percentage of Interaction (%)'].max()
        data['frequency'] = (data['Percentage of Interaction (%)'] - (min_percentage - 1)) / ((1 + max_percentage) - (min_percentage - 1))
        min_energy = data['Lowest Energy'].min()
        max_energy = data['Lowest Energy'].max()
        data['Energy'] = ((1 + max_energy) - data['Lowest Energy']) / ((1 + max_energy) - (min_energy - 1))
        data['FP_score'] = data['frequency'] * data['Energy']
        data.to_csv('FP_scores.csv', index=False)
        file_path = 'FP_scores.csv'
        df = pd.read_csv(file_path, sep=',')
        df['Formatted_Fragment'] = df['Fragment'].apply(format_fragment_name)
        unique_fragments = df['Formatted_Fragment'].unique()
        palette = sns.color_palette("hsv", len(unique_fragments))
        color_map = dict(zip(unique_fragments, palette))
        mean_frequency = df['frequency'].mean()
        mean_energy = df['Energy'].mean()
        above_threshold_points = df[(df['frequency'] > mean_frequency) | (df['Energy'] > mean_energy)]
        if not above_threshold_points.empty:
            threshold = above_threshold_points['FP_score'].min()
        else:
            threshold = 0.1

        fig1, ax1 = plt.subplots(figsize=(8, 5))  # Reduced size but keep font size high

        ax1.fill_betweenx([0, mean_energy], mean_frequency, 1, color='green', alpha=0.1)  
        ax1.fill_betweenx([mean_energy, 1], 0, mean_frequency, color='yellow', alpha=0.1)  
        ax1.fill_betweenx([mean_energy, 1], mean_frequency, 1, color='red', alpha=0.1)  
        ax1.fill_betweenx([0, mean_energy], 0, mean_frequency, color='cyan', alpha=0.1) 

        texts = []
        for fragment in unique_fragments:
            subset = df[df['Formatted_Fragment'] == fragment]
            ax1.scatter(subset['frequency'], subset['Energy'], color=color_map[fragment], s=100, edgecolor='black') 
            for i, row in subset.iterrows():
                texts.append(ax1.text(row['frequency'], row['Energy'], row['Formatted_Fragment'], fontsize=14))  # Fontsize reduced for text
        adjust_text(texts, ax=ax1, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))
        # ax1.axhline(mean_energy, color='grey', linestyle='--', label=f'Mean Energy: {mean_energy:.2f}')
        # ax1.axvline(mean_frequency, color='grey', linestyle='--', label=f'Mean Frequency: {mean_frequency:.2f}')
        ax1.axhline(mean_energy, color='grey', linestyle='--')
        ax1.axvline(mean_frequency, color='grey', linestyle='--')
        ax1.set_xlabel('Normalization of Percentage of interaction (PI%)', fontsize=14, fontweight='bold', color='black')  # Reduced fontsize
        ax1.set_ylabel('Normalization of Interaction Energy (ΔE)', fontsize=14, fontweight='bold', color='black')  # Reduced fontsize
        # ax1.set_title('Energy vs Percentage of Interaction', fontsize=16, fontweight='bold')  # Slightly larger for the title
        # ax1.legend(title='Mean Values')
        ax1.tick_params(axis='x', labelsize=12)  # Reduced x-axis tick label font size
        ax1.tick_params(axis='y', labelsize=12)
        plt.setp(ax1.get_xticklabels(), fontweight='bold')
        plt.setp(ax1.get_yticklabels(), fontweight='bold')
        plt.tight_layout()
        plt.savefig('library_analysis/FP_score_vs_Hotspots.png', dpi=300)  # High DPI to maintain quality at smaller size
        fig2, ax2 = plt.subplots(figsize=(8, 5))  # Reduced size but keep font size high
        colors = [color_map[fragment] for fragment in df['Formatted_Fragment']]
        ax2.bar(df['Formatted_Fragment'], df['FP_score'], width=0.4, alpha=0.6, color=colors, edgecolor='black')
        ax2.axhline(threshold, color='black', linestyle='--', label=f'Threshold: {threshold:.2f}')
        ax2.set_xlabel('Binding site Hotspots', fontsize=14, fontweight='bold')
        ax2.set_ylabel('FP_score', fontsize=14, fontweight='bold', color='black')
        ax2.tick_params(axis='y', labelcolor='black', labelsize=12)
        ax2.tick_params(axis='x', labelsize=12, rotation=90)
        ax2.set_title('FP-score vs Binding site Hotspots', fontsize=16, fontweight='bold')
        ax2.legend()
        plt.setp(ax2.get_xticklabels(), fontweight='bold')
        plt.setp(ax2.get_yticklabels(), fontweight='bold')
        plt.tight_layout()
        plt.savefig('library_analysis/Energy_vs_Percentage_of_Interaction.png', dpi=300)
        ##########################################################################################################################
        file_path = 'library_analysis/merged_file.csv'  
        df = pd.read_csv(file_path)
        df['Formatted_Fragment'] = df['Fragment'].apply(format_fragment_name)
        # Define bond types
        bond_types = [
            "Halogen_Bonds", 
            "Hydrogen_Bonds", 
            "Hydrophobic_Interactions", 
            "Water_Bridges", 
            "pi-Stacking",
            "pi-Cation_Interactions",
            "Metal_Complexes",
            "Salt_Bridges"
        ]

        for bond_type in bond_types:
            if bond_type not in df.columns:
                df[bond_type] = 0

        # Normalize bond types to match "Percentage of Interaction" for each fragment
        for i, row in df.iterrows():
            total_bond_interactions = row[bond_types].sum()
            if total_bond_interactions > 0:
                df.loc[i, bond_types] = (row[bond_types] / total_bond_interactions) * row["Percentage of Interaction (%)"]

        # Create the plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [1, 1]})
        bottom = np.zeros(len(df))

        color_mapping = {
            'pi-Stacking': 'orange',
            'Hydrophobic_Interactions': 'green',
            'Hydrogen_Bonds': 'blue',
            'Halogen_Bonds': 'cyan',
            'Salt_Bridges': 'red',
            'Water_Bridges': 'magenta',
            'pi-Cation_Interactions': 'yellow',
            'Metal_Complexes': 'brown'
        }

        # Stack the bars for different bond types
        for bond_type in bond_types:
            color = color_mapping.get(bond_type, 'gray')
            ax1.bar(df["Formatted_Fragment"], df[bond_type], bottom=bottom, label=bond_type, color=color, edgecolor='black')
            bottom += df[bond_type]

        # ax1.set_ylabel("Percentage of Interaction (%)", fontsize=16, fontweight='bold')  # Increased font size
        ax1.set_ylim(0, 100)
        ax1.tick_params(axis='y', labelsize=14, labelcolor='black')
        # ax1.set_title("Percentage of Interaction by Bond Type", fontsize=18, fontweight='bold')  # Add title
        # ax1.legend(fontsize=12, title='Bond Type', title_fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        norm = plt.Normalize(-50, 20)
        colors = plt.cm.coolwarm(norm(df["Lowest Energy"]))

        ax2.bar(df["Formatted_Fragment"], df["Lowest Energy"], color=colors, edgecolor='black')
        # ax2.set_ylabel("Lowest Energy", fontsize=16, fontweight='bold')  # Increased font size
        ax2.set_ylim(20, -50)
        ax2.invert_yaxis()
        ax2.tick_params(axis='y', labelsize=14, labelcolor='black')
        sm = plt.cm.ScalarMappable(cmap=cm.coolwarm, norm=norm)
        sm.set_array([]) 
        cbar = plt.colorbar(sm, ax=ax2, orientation='horizontal', pad=0.4)
        # cbar.set_label('Lowest Energy', fontsize=16, fontweight='bold')  # Increased font size
        cbar.ax.tick_params(labelsize=12)
        plt.xticks(rotation=90, fontsize=14, fontweight='bold')  # Increased x-axis tick label font size
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust rect to make space for the legend
        output_file_coolwarm = 'library_analysis/Binding_site_residues_barchart_corrected.png'
        plt.savefig(output_file_coolwarm, dpi=400) 
        os.remove('library_analysis/binding_site_residues_bond_type_percentages.csv')
        os.remove('library_analysis/binding_site_residues_interaction_and_energy.csv')
        # ###############################################################################
        parent_dir = os.getcwd()
        all_LLEs = []

        def LLE_def(parent_dir):
            all_LLEs = []
            for subdir in os.listdir(parent_dir):
                subdir_path = os.path.join(parent_dir, subdir)
                if not os.path.isdir(subdir_path):
                    continue
                for subdir_chain in os.listdir(subdir_path):
                    subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                    if not os.path.isdir(subdir_chain_path):
                        continue
                    analysis_dir = os.path.join(subdir_chain_path, "analysis")
                    if not os.path.isdir(analysis_dir):
                        continue
                    os.chdir(analysis_dir)
                    if not os.path.isfile("LLE.txt"):
                        continue
                    with open('LLE.txt', 'r') as infile:
                        for line in infile:
                            if line.startswith("total_LLE"):
                                LLE = float(line.split('=')[1].strip())
                                all_LLEs.append((subdir_chain, round(LLE, 1)))
            os.chdir(parent_dir)

            sorted_LLEs = sorted(all_LLEs, key=lambda x: x[1], reverse=True)
            output_dir = os.path.join(parent_dir, 'library_analysis')
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, 'LLE.txt')
            with open(output_file, 'w') as outfile:
                outfile.write(f"Fragment   LLE\n")
                for folder, LLE in sorted_LLEs:
                    outfile.write(f"{folder}    {LLE}\n")

            print(f"Data has been processed and saved to {output_file}")

        parent_dir = os.getcwd() 
        LLE_def(parent_dir)

        ###################################################################################################
        def process_fragments(parent_dir):
            output_file = os.path.join("library_analysis", "LLE_all.txt")
            suggested_frags = os.path.join(parent_dir, "library_analysis/suggested_frags")
            os.makedirs(suggested_frags, exist_ok=True)
            all_fragments = []

            for subdir in os.listdir(parent_dir):
                subdir_path = os.path.join(parent_dir, subdir)
                if not os.path.isdir(subdir_path):
                    continue
                for subdir_chain in os.listdir(subdir_path):
                    subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                    if not os.path.isdir(subdir_chain_path):
                        continue
                    analysis_dir = os.path.join(subdir_chain_path, "analysis")
                    if not os.path.isdir(analysis_dir):
                        continue
                    lle_file_path = os.path.join(analysis_dir, "LLE.txt")
                    if not os.path.isfile(lle_file_path):
                        continue
                    with open(lle_file_path, 'r') as file:
                        lines = file.readlines()
                        for line in lines[1:]:
                            line = line.strip()
                            if not line or line.startswith("total_"):
                                continue
                            parts = line.split(maxsplit=4)
                            if len(parts) < 5:
                                continue
                            fragment = parts[0]
                            if fragment.startswith("linker") or fragment.startswith("RingSystem"):
                                try:
                                    total_energy = float(parts[1])
                                    lle = float(parts[4].split()[0])
                                    hotspot = parts[4][parts[4].find("["):].strip() if "[" in parts[4] else ""
                                    all_fragments.append([f"{subdir}/{subdir_chain}/{fragment}", total_energy, lle, hotspot])
                                except (ValueError, IndexError) as e:
                                    print(f"Error processing line: {line}, error: {e}")

            df = pd.DataFrame(all_fragments, columns=["Fragment", "total_energy", "LLE", "Hotspot"])
            df_sorted = df.sort_values(by="LLE", ascending=False)

            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as file:
                file.write("Fragment total_energy LLE Hotspot\n")
                for index, row in df_sorted.iterrows():
                    file.write(f"{row['Fragment']} {row['total_energy']} {row['LLE']} {row['Hotspot']}\n")

            print(f"LLE_all has been processed and saved to {output_file}")

            unique_residues = {}
            for index, row in df.iterrows():
                residue = row["Hotspot"]
                if residue not in unique_residues:
                    unique_residues[residue] = row
                elif row["total_energy"] < unique_residues[residue]["total_energy"]:
                    unique_residues[residue] = row

            df_unique = pd.DataFrame(unique_residues.values(), columns=df.columns)
            for index, row in df_unique.iterrows():
                fragment_path = row['Fragment']
                fragment_dir, fragment_name = os.path.split(fragment_path)
                subdir_chain = fragment_dir.split('/')[-1]
                pdb_file = f"{fragment_name}.pdb"
                source_path = os.path.join(parent_dir, fragment_dir, "analysis", pdb_file)
                if os.path.isfile(source_path):
                    dest_file = f"{subdir_chain.replace('/', '_')}_{fragment_name}.pdb"
                    dest_path = os.path.join(suggested_frags, dest_file)
                    shutil.copyfile(source_path, dest_path)

            txt_file = os.path.join(parent_dir, "library_analysis/Ph4_lig.pdb")
            matching_dir = os.path.join(parent_dir, "library_analysis/matching_frags")
            os.makedirs(matching_dir, exist_ok=True)

            def extract_coordinates(line):
                return line[30:54]

            coordinates_set = set()
            with open(txt_file, 'r') as file:
                for line in file:
                    coordinates_set.add(extract_coordinates(line))

            def has_matching_coordinates(pdb_file, coordinates_set):
                with open(pdb_file, 'r') as file:
                    for line in file:
                        if extract_coordinates(line) in coordinates_set:
                            return True
                return False

            for pdb_filename in os.listdir(suggested_frags):
                if pdb_filename.endswith(".pdb"):
                    pdb_path = os.path.join(suggested_frags, pdb_filename)
                    if has_matching_coordinates(pdb_path, coordinates_set):
                        matching_path = os.path.join(matching_dir, pdb_filename)
                        shutil.copyfile(pdb_path, matching_path)
                        print(f"Matching coordinates found in {pdb_filename}")

        parent_dir = os.getcwd()
        print(parent_dir)
        process_fragments(parent_dir)
        ###################################################################################################
        data = pd.read_csv('library_analysis/merged_file.csv')
        min_percentage = data['Percentage of Interaction (%)'].min()
        max_percentage = data['Percentage of Interaction (%)'].max()
        data['frequency'] = (data['Percentage of Interaction (%)'] - (min_percentage - 1)) / ((1 + max_percentage) - (min_percentage - 1))
        min_energy = data['Lowest Energy'].min()
        max_energy = data['Lowest Energy'].max()
        data['Energy'] = ((1 + max_energy) - data['Lowest Energy']) / ((1 + max_energy) - (min_energy - 1))
        data['FP_score'] = data['frequency'] * data['Energy']
        data.to_csv('library_analysis/FP_scores.csv', index=False)
        file_path = 'library_analysis/FP_scores.csv'
        df = pd.read_csv(file_path, sep=',')
        df['Formatted_Fragment'] = df['Fragment'].apply(format_fragment_name)
        unique_fragments = df['Formatted_Fragment'].unique()
        palette = sns.color_palette("hsv", len(unique_fragments))
        color_map = dict(zip(unique_fragments, palette))
        mean_frequency = df['frequency'].mean()
        mean_energy = df['Energy'].mean()
        above_threshold_points = df[(df['frequency'] > mean_frequency) | (df['Energy'] > mean_energy)]
        if not above_threshold_points.empty:
            threshold = above_threshold_points['FP_score'].min()
        else:
            threshold = 0.1
        fig1, ax1 = plt.subplots(figsize=(15, 8))
        ax1.fill_betweenx([0, mean_energy], mean_frequency, 1, color='green', alpha=0.1)  
        ax1.fill_betweenx([mean_energy, 1], 0, mean_frequency, color='yellow', alpha=0.1)  
        ax1.fill_betweenx([mean_energy, 1], mean_frequency, 1, color='red', alpha=0.1)  
        ax1.fill_betweenx([0, mean_energy], 0, mean_frequency, color='cyan', alpha=0.1) 
        texts = []
        for fragment in unique_fragments:
            subset = df[df['Formatted_Fragment'] == fragment]
            ax1.scatter(subset['frequency'], subset['Energy'], color=color_map[fragment], s=100, edgecolor='black') 
            for i, row in subset.iterrows():
                texts.append(ax1.text(row['frequency'], row['Energy'], row['Formatted_Fragment'], fontsize=20))  # Use formatted fragment names
        adjust_text(texts, ax=ax1, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))

        # Add mean lines
        ax1.axhline(mean_energy, color='grey', linestyle='--', label=f'Mean Energy: {mean_energy:.2f}')
        ax1.axvline(mean_frequency, color='grey', linestyle='--', label=f'Mean Frequency: {mean_frequency:.2f}')

        # Set labels and titles for the scatter plot
        ax1.set_xlabel('Normalization of Percentage of interaction (PI%)', fontsize=20, fontweight='bold', color='black')
        ax1.set_ylabel('Normalization of Interaction Energy (ΔE)', fontsize=20, fontweight='bold', color='black')
        ax1.set_title('Energy vs Percentage of Interaction', fontsize=20, fontweight='bold')
        ax1.legend(title='Mean Values')

        # Bold the tick labels
        ax1.tick_params(axis='x', labelsize=18)  # Increase x-axis tick label font size
        ax1.tick_params(axis='y', labelsize=18)

        plt.setp(ax1.get_xticklabels(), fontweight='bold')
        plt.setp(ax1.get_yticklabels(), fontweight='bold')

        # Save the first plot to a PNG file
        plt.tight_layout()
        plt.savefig('library_analysis/FP_score_vs_Hotspots.png')
        # Plot 2: P-score vs Binding site Hotspots (bar chart)
        fig2, ax2 = plt.subplots(figsize=(15, 8))

        # Plot the P-scores on the bar chart
        colors = [color_map[fragment] for fragment in df['Formatted_Fragment']]
        ax2.bar(df['Formatted_Fragment'], df['FP_score'], width=0.4, alpha=0.6, color=colors, edgecolor='black')

        # Add the threshold line
        ax2.axhline(threshold, color='black', linestyle='--', label=f'Threshold: {threshold:.2f}')

        # Set labels and titles for the bar chart
        ax2.set_xlabel('Binding site Hotspots', fontsize=20, fontweight='bold')
        ax2.set_ylabel('FP_score', fontsize=20, fontweight='bold', color='black')
        ax2.tick_params(axis='y', labelcolor='black', labelsize=20)
        ax2.tick_params(axis='x', labelsize=20, rotation=90)
        ax2.set_title('FP-score vs Binding site Hotspots', fontsize=20, fontweight='bold')
        ax2.legend()

        # Bold the tick labels
        plt.setp(ax2.get_xticklabels(), fontweight='bold')
        plt.setp(ax2.get_yticklabels(), fontweight='bold')
        plt.tight_layout()
        plt.savefig('library_analysis/Energy_vs_Percentage_of_Interaction.png')
        ###################################################################################################


        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        from adjustText import adjust_text
        import numpy as np
        import seaborn as sns
        import matplotlib.cm as cm
        import os

        def format_fragment_name(fragment):
            residue = fragment[:3].capitalize()
            chain = fragment[3]
            number = fragment[4:]
            return f"{residue}-{chain}-{number}"
        data = pd.read_csv('library_analysis/merged_file.csv')
        min_percentage = data['Percentage of Interaction (%)'].min()
        max_percentage = data['Percentage of Interaction (%)'].max()
        data['frequency'] = (data['Percentage of Interaction (%)'] - (min_percentage - 1)) / ((1 + max_percentage) - (min_percentage - 1))
        min_energy = data['Lowest Energy'].min()
        max_energy = data['Lowest Energy'].max()
        data['Energy'] = ((1 + max_energy) - data['Lowest Energy']) / ((1 + max_energy) - (min_energy - 1))
        data['FP_score'] = data['frequency'] * data['Energy']
        data.to_csv('library_analysis/FP_scores.csv', index=False)
        file_path = 'library_analysis/FP_scores.csv'
        df = pd.read_csv(file_path, sep=',')
        df['Formatted_Fragment'] = df['Fragment'].apply(format_fragment_name)
        unique_fragments = df['Formatted_Fragment'].unique()
        palette = sns.color_palette("hsv", len(unique_fragments))
        color_map = dict(zip(unique_fragments, palette))
        mean_frequency = df['frequency'].mean()
        mean_energy = df['Energy'].mean()
        above_threshold_points = df[(df['frequency'] > mean_frequency) | (df['Energy'] > mean_energy)]
        if not above_threshold_points.empty:
            threshold = above_threshold_points['FP_score'].min()
        else:
            threshold = 0.1
        fig1, ax1 = plt.subplots(figsize=(8, 5))
        ax1.fill_betweenx([0, mean_energy], mean_frequency, 1, color='green', alpha=0.1)  
        ax1.fill_betweenx([mean_energy, 1], 0, mean_frequency, color='yellow', alpha=0.1)  
        ax1.fill_betweenx([mean_energy, 1], mean_frequency, 1, color='red', alpha=0.1)  
        ax1.fill_betweenx([0, mean_energy], 0, mean_frequency, color='cyan', alpha=0.1) 
        texts = []
        for fragment in unique_fragments:
            subset = df[df['Formatted_Fragment'] == fragment]
            ax1.scatter(subset['frequency'], subset['Energy'], color=color_map[fragment], s=100, edgecolor='black') 
            for i, row in subset.iterrows():
                texts.append(ax1.text(row['frequency'], row['Energy'], row['Formatted_Fragment'], fontsize=14))  
        adjust_text(texts, ax=ax1, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))
        ax1.axhline(mean_energy, color='grey', linestyle='--')
        ax1.axvline(mean_frequency, color='grey', linestyle='--')
        ax1.set_xlabel('Normalization of Percentage of interaction (PI%)', fontsize=14, fontweight='bold', color='black')  
        ax1.set_ylabel('Normalization of Interaction Energy (ΔE)', fontsize=14, fontweight='bold', color='black')  
        ax1.tick_params(axis='x', labelsize=12)
        ax1.tick_params(axis='y', labelsize=12)
        plt.setp(ax1.get_xticklabels(), fontweight='bold')
        plt.setp(ax1.get_yticklabels(), fontweight='bold')
        plt.tight_layout()
        plt.savefig('library_analysis/FP_score_vs_Hotspots.png', dpi=300)  
        fig2, ax2 = plt.subplots(figsize=(8, 5))  
        colors = [color_map[fragment] for fragment in df['Formatted_Fragment']]
        ax2.bar(df['Formatted_Fragment'], df['FP_score'], width=0.4, alpha=0.6, color=colors, edgecolor='black')
        ax2.axhline(threshold, color='black', linestyle='--', label=f'Threshold: {threshold:.2f}')
        ax2.set_xlabel('Binding site Hotspots', fontsize=14, fontweight='bold')
        ax2.set_ylabel('FP_score', fontsize=14, fontweight='bold', color='black')
        ax2.tick_params(axis='y', labelcolor='black', labelsize=12)
        ax2.tick_params(axis='x', labelsize=12, rotation=90)
        ax2.set_title('FP-score vs Binding site Hotspots', fontsize=16, fontweight='bold')
        ax2.legend()
        plt.setp(ax2.get_xticklabels(), fontweight='bold')
        plt.setp(ax2.get_yticklabels(), fontweight='bold')
        plt.tight_layout()
        plt.savefig('library_analysis/Energy_vs_Percentage_of_Interaction.png', dpi=300) 
        file_path = 'library_analysis/merged_file.csv'  
        df = pd.read_csv(file_path)
        df['Formatted_Fragment'] = df['Fragment'].apply(format_fragment_name)
        bond_types = [
            "Halogen_Bonds", 
            "Hydrogen_Bonds", 
            "Hydrophobic_Interactions", 
            "Water_Bridges", 
            "pi-Stacking",
            "pi-Cation_Interactions",
            "Metal_Complexes",
            "Salt_Bridges"
        ]
        for bond in bond_types:
            if bond not in df.columns:
                df[bond] = 0
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [1, 1]})
        bottom = np.zeros(len(df))
        color_mapping = {
            'pi-Stacking': 'orange',
            'Hydrophobic_Interactions': 'green',
            'Hydrogen_Bonds': 'blue',
            'Halogen_Bonds': 'cyan',
            'Salt_Bridges': 'red',
            'Water_Bridges': 'magenta',
            'pi-Cation_Interactions': 'yellow',
            'Metal_Complexes': 'brown'
        }
        for bond_type in bond_types:
            color = color_mapping.get(bond_type, 'gray')
            ax1.bar(df["Formatted_Fragment"], df[bond_type], bottom=bottom, label=bond_type, color=color, edgecolor='black')
            bottom += df[bond_type]
        ax1.set_ylim(0, 30)
        ax1.tick_params(axis='y', labelsize=14, labelcolor='black')
        norm = plt.Normalize(-50, 20)
        colors = plt.cm.coolwarm(norm(df["Lowest Energy"]))
        ax2.bar(df["Formatted_Fragment"], df["Lowest Energy"], color=colors, edgecolor='black')
        ax2.set_ylim(20, -50)
        ax2.invert_yaxis()
        ax2.tick_params(axis='y', labelsize=14, labelcolor='black')
        sm = plt.cm.ScalarMappable(cmap=cm.coolwarm, norm=norm)
        sm.set_array([]) 
        cbar = plt.colorbar(sm, ax=ax2, orientation='horizontal', pad=0.4)
        cbar.ax.tick_params(labelsize=12)
        plt.xticks(rotation=90, fontsize=14, fontweight='bold')  
        plt.tight_layout(rect=[0, 0, 0.85, 1])  
        output_file_coolwarm = 'library_analysis/Binding_site_residues_barchart_corrected.png'
        plt.savefig(output_file_coolwarm, dpi=400) 
        parent_dir = os.getcwd()
        selected_data_path = os.path.join(parent_dir, 'library_analysis/selected_data.csv')
        selected_data = pd.read_csv(selected_data_path, index_col=0)
        selected_data['total'] = selected_data.sum(axis=1)
        selected_data = selected_data.sort_values('total', ascending=True)
        selected_data = selected_data.drop('total', axis=1)
        def get_chain_id_and_residue_number(column_name):
            chain_id = column_name[3]
            residue_number = column_name[4:]
            return (chain_id, int(residue_number))
        sorted_columns = sorted(selected_data.columns, key=get_chain_id_and_residue_number)
        selected_data = selected_data[sorted_columns]
        percentage_found = (selected_data.count() / len(selected_data) * 100).round(0).astype(int)
        formatted_columns = [format_fragment_name(col) for col in sorted_columns]
        plt.style.use('seaborn-darkgrid')
        plt.figure(figsize=(12, 6))  
        ax = sns.heatmap(selected_data,
                         linewidths=0.5,
                         linecolor='grey',
                         cmap='coolwarm',
                         center=0,
                         annot=True,
                         annot_kws={"size": 8},
                         fmt=".1f",
                         vmax=30,
                         vmin=-50)
        ax.vlines(range(len(selected_data.columns) + 1), *ax.get_ylim(), colors='grey', linewidth=0.5)
        ax.set_facecolor('white')
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        tick_positions = np.arange(0.5, len(selected_data.columns), 1)
        ax2.set_xticks(tick_positions)
        ax2.set_xticklabels(percentage_found, ha='center', rotation=0, fontsize=12, fontweight='bold')
        ax2.grid(False) 
        plt.tick_params(labelsize=12)
        plt.yticks(rotation='horizontal', fontweight='bold')  
        plt.xticks(rotation=0, ha="center", fontweight='bold')  
        plt.setp(ax.get_xticklabels(), fontsize=12, fontweight='bold')
        plt.setp(ax.get_yticklabels(), fontsize=12, fontweight='bold')
        ax.set_xticklabels(formatted_columns)
        heatmap_filename = os.path.join(parent_dir, 'library_analysis/Ph4_heatmap.png')
        plt.tight_layout()
        plt.savefig(heatmap_filename, dpi=400)

        print(os.getcwd())
        # output_file = "../../../library_analysis/pocket.pdb"
        def process_pocket_files(parent_dir):
            output_file = "library_analysis/pocket.pdb"
            written_lines = set()

            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, "w") as output_f:
                for subdir in os.listdir(parent_dir):
                    subdir_path = os.path.join(parent_dir, subdir)
                    if not os.path.isdir(subdir_path):
                        continue
                    for subdir_chain in os.listdir(subdir_path):
                        subdir_chain_path = os.path.join(subdir_path, subdir_chain)
                        if not os.path.isdir(subdir_chain_path):
                            continue
                        analysis_BS_path = os.path.join(subdir_chain_path, "FMOPhore")
                        if not os.path.isdir(analysis_BS_path):
                            continue
                        for filename in os.listdir(analysis_BS_path):
                            if filename != "pocket.pdb":
                                continue
                            file_path = os.path.join(analysis_BS_path, filename)
                            with open(file_path, "r") as input_f:
                                lines1 = input_f.readlines()
                                lines1 = [line for line in lines1 if 'HEADER' not in line]
                                lines1 = [line for line in lines1 if 'END' not in line]
                                lines1 = [line for line in lines1 if 'TER' not in line]
                                for line in lines1:
                                    key = line[11:26]
                                    if key not in written_lines:
                                        output_f.write(line)
                                        written_lines.add(key)

            print(f"Data from all pocket.pdb files has been combined and saved to {output_file}")
        parent_dir = os.getcwd() 
        process_pocket_files(parent_dir)

        import pandas as pd
        from pymol import cmd
        import os
        import glob
        import shutil 

        parent_folder = './'  
        target_folder = os.path.join(parent_folder, 'library_pymol')
        os.makedirs(target_folder, exist_ok=True)
        csv_file = os.path.join(parent_folder, 'FP_scores.csv')
        FP_scores = pd.read_csv(csv_file)
        mean_energy = FP_scores['Energy'].mean()
        mean_frequency = FP_scores['frequency'].mean()
        def process_fmophore(fmophore_folder):
            ph4_3d_fmo_file = os.path.join(fmophore_folder, 'Ph4_3D_FMO.txt')
            if not os.path.exists(ph4_3d_fmo_file):
                print(f"File {ph4_3d_fmo_file} not found. Skipping...")
                return
            Ph4_3D_FMO = pd.read_csv(ph4_3d_fmo_file, sep=r'[\t\s]+', engine='python')
            Ph4_3D_FMO.columns = Ph4_3D_FMO.columns.str.strip()
            Ph4_3D_FMO = Ph4_3D_FMO.applymap(lambda x: x.strip() if isinstance(x, str) else x)
            for col in ['LIG.x', 'LIG.y', 'LIG.z', 'PROT.x', 'PROT.y', 'PROT.z']:
                Ph4_3D_FMO[col] = pd.to_numeric(Ph4_3D_FMO[col], errors='coerce')
            color_mapping = {
                'pi-Stacking': 'orange',
                'Hydrophobic_Interactions': 'green',
                'Hydrogen_Bonds': 'blue',
                'Halogen_Bonds': 'cyan',
                'Salt_Bridges': 'red',
                'Water_Bridges': 'magenta',
                'pi-Cation_Interactions': 'yellow',
                'Metal_Complexes': 'brown'
            }
            pdb_file_pattern = os.path.join(fmophore_folder, '../*_FMO_protein.pdb')
            pdb_lig_file_pattern = os.path.join(fmophore_folder, '../*_FMO_lig_H.pdb')
            ph4_lig_file = os.path.join(fmophore_folder, 'Ph4_lig.pdb')  # Added Ph4_lig.pdb
            pdb_file = glob.glob(pdb_file_pattern)
            pdb_lig_file = glob.glob(pdb_lig_file_pattern)
            print(f"PDB files found: {pdb_file}")
            print(f"Ligand files found: {pdb_lig_file}")
            if not pdb_file or not pdb_lig_file:
                print(f"PDB or ligand file not found in {fmophore_folder}. Skipping...")
                return
            pdb_file = pdb_file[0]
            pdb_lig_file = pdb_lig_file[0]
            pdb_base_name = os.path.splitext(os.path.basename(pdb_file))[0]
            lig_base_name = os.path.splitext(os.path.basename(pdb_lig_file))[0]
            cmd.load(pdb_lig_file, f'{lig_base_name}_ligand')
            cmd.load(pdb_file, f'{pdb_base_name}_protein')
            cmd.bg_color('white')
            cmd.show("cartoon", f'{pdb_base_name}_protein')
            cmd.color("white", f'{pdb_base_name}_protein')
            cmd.hide("everything", "(hydro and elem H and neighbor elem C)")
            cmd.set('cartoon_transparency', 0.8)
            group_name = f'{pdb_base_name}_ligand_interactions_group'
            cmd.color('magenta', f'{lig_base_name}_ligand')
            cmd.group(group_name, f'{lig_base_name}_ligand')
            cmd.group(group_name, f'{pdb_base_name}_protein')
            for index, row in Ph4_3D_FMO.iterrows():
                fragment = row['Fragment']
                chain_id = fragment[3]
                residue_number = fragment[4:]
                selection_name = f"{residue_number}{chain_id}"
                selection_command = f"resi {residue_number} and chain {chain_id}"
                lig_x, lig_y, lig_z = row['LIG.x'], row['LIG.y'], row['LIG.z']
                prot_x, prot_y, prot_z = row['PROT.x'], row['PROT.y'], row['PROT.z']
                if pd.notna([lig_x, lig_y, lig_z, prot_x, prot_y, prot_z]).all():
                    lig_sphere_name = f'{pdb_base_name}_lig_sphere_{index}'
                    prot_sphere_name = f'{pdb_base_name}_prot_sphere_{index}'
                    cmd.pseudoatom(lig_sphere_name, pos=[lig_x, lig_y, lig_z], color="magenta")
                    cmd.pseudoatom(prot_sphere_name, pos=[prot_x, prot_y, prot_z], color="yellow")
                    cmd.group(group_name, lig_sphere_name)
                    cmd.group(group_name, prot_sphere_name)
                    dash_name = f'{pdb_base_name}_dash_{index}'
                    cmd.distance(dash_name, lig_sphere_name, prot_sphere_name)
                    cmd.hide('labels', f'{dash_name}')
                    cmd.group(group_name, dash_name)
                    cmd.set('dash_radius', 0.05)
                    bond_type = row['Bond_type'].strip()
                    if bond_type in color_mapping:
                        dash_color = color_mapping[bond_type]
                        cmd.color(dash_color, dash_name)
                    else:
                        print(f"Warning: Bond type '{bond_type}' not found in color mapping.")
                matching_row = FP_scores[FP_scores['Fragment'] == fragment]
                if not matching_row.empty:
                    energy = matching_row['Energy'].values[0]
                    frequency = matching_row['frequency'].values[0]
                    cmd.select(selection_name, selection_command)
                    cmd.show("sticks", selection_name)
                    cmd.hide('sticks', f'{selection_name} and elem H and neighbor elem C')
                    if energy < mean_energy and frequency > mean_frequency:
                        cmd.color("limegreen", selection_name)
                    elif energy > mean_energy and frequency < mean_frequency:
                        cmd.color("yelloworange", selection_name)
                    elif energy > mean_energy and frequency > mean_frequency:
                        cmd.color("deepsalmon", selection_name)
                    elif energy < mean_energy and frequency < mean_frequency:
                        cmd.color("aquamarine", selection_name)
                    cmd.group(group_name, selection_name)
            if os.path.exists(ph4_lig_file):
                ph4_lig_name = os.path.splitext(os.path.basename(ph4_lig_file))[0]
                cmd.load(ph4_lig_file, f'{ph4_lig_name}')
                cmd.hide('sticks', f'{ph4_lig_name}') 
                cmd.show('spheres', f'{ph4_lig_name}')
                cmd.set('sphere_scale', 0.3, f'{ph4_lig_name}')  
                cmd.set('sphere_transparency', 0.3, f'{ph4_lig_name}')  
                cmd.color('green', f'{ph4_lig_name} and elem C')  
                cmd.color('white', f'{ph4_lig_name} and elem H')  
                cmd.color('blue', f'{ph4_lig_name} and elem N')  
                cmd.color('red', f'{ph4_lig_name} and elem O')    
                cmd.color('yellow', f'{ph4_lig_name} and elem S') 
                cmd.group(group_name, f'{ph4_lig_name}') 
            pymol_session_path = os.path.join(fmophore_folder, f'{pdb_base_name}_pymol.pse')
            cmd.zoom()
            cmd.save(pymol_session_path)
            print(f"PyMOL session saved as {pymol_session_path}!")
            target_session_path = os.path.join(target_folder, f'{pdb_base_name}_pymol.pse')
            shutil.move(pymol_session_path, target_session_path)
            print(f"Moved {pdb_base_name}_pymol_session_with_spheres.pse to {target_session_path}!")
            cmd.reinitialize()
        def process_suggested_hotspots(fmophore_folder):
            pdb_file_pattern = os.path.join(fmophore_folder, '../*_FMO_protein.pdb')
            pdb_files = glob.glob(pdb_file_pattern)
            if pdb_files:
                pdb_file = pdb_files[0]  
                pdb_base_name = os.path.basename(pdb_file).split('.')[0]
                cmd.load(pdb_file, f'{pdb_base_name}_protein')
                cmd.bg_color('white')
                cmd.show("cartoon", f'{pdb_base_name}_protein')
                cmd.color("white", f'{pdb_base_name}_protein')
                cmd.hide("everything", "(hydro and elem H and neighbor elem C)")
                cmd.set('cartoon_transparency', 0.8)
                group_name = f"{pdb_base_name}_group"
                cmd.group(group_name, f'{pdb_base_name}_protein')
                for index, row in FP_scores.iterrows():
                    fragment = row['Fragment']
                    chain_id = fragment[3]
                    residue_number = fragment[4:]
                    selection_name = f"{residue_number}{chain_id}"
                    selection_command = f"resi {residue_number} and chain {chain_id}"
                    matching_row = FP_scores[FP_scores['Fragment'] == fragment]
                    if not matching_row.empty:
                        energy = matching_row['Energy'].values[0]
                        frequency = matching_row['frequency'].values[0]
                        cmd.select(selection_name, selection_command)
                        cmd.show("cartoon", selection_name)  
                        cmd.show("sticks", selection_name)  
                        cmd.hide('sticks', f'{selection_name} and elem H and neighbor elem C')
                        if energy < mean_energy and frequency > mean_frequency:
                            cmd.color("limegreen", selection_name)
                        elif energy > mean_energy and frequency < mean_frequency:
                            cmd.color("yelloworange", selection_name)
                        elif energy > mean_energy and frequency > mean_frequency:
                            cmd.color("deepsalmon", selection_name)
                        elif energy < mean_energy and frequency < mean_frequency:
                            cmd.color("aquamarine", selection_name)
                        cmd.group(group_name, selection_name)
                pymol_session_path = os.path.join(fmophore_folder, f'{pdb_base_name}_suggested_pymol.pse')
                cmd.zoom()
                cmd.save(pymol_session_path)
                print(f"PyMOL session saved as {pymol_session_path}!")
                target_session_path = os.path.join(target_folder, f'{pdb_base_name}_suggested_pymol.pse')
                shutil.move(pymol_session_path, target_session_path)
                print(f"Moved {pdb_base_name}_pymol.pse to {target_session_path}!")
                cmd.reinitialize()
        fmophore_folders = glob.glob(os.path.join(parent_folder, '**/FMOPhore'), recursive=True)
        for folder in fmophore_folders:
            print(f"Processing folder: {folder}")
            process_fmophore(folder)
            process_suggested_hotspots(folder)
        ###############################################################################
        def safe_remove(file_path):
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    print(f"Removed: {file_path}")
            except Exception as e:
                print(f"Error removing {file_path}: {e}")
        def main(parent_dir):
            print(f"Current working directory: {os.getcwd()}")
            library_analysis_path = os.path.join(parent_dir, "library_analysis")
            common_files_to_remove = ["pocket.pqr", "time.txt", "mdpout_dens_grid.dx", "mdpout_freq_grid.dx", "mdpout_all_atom_pdensities.pdb"]
            for file_name in os.listdir(library_analysis_path):
                if file_name.endswith('.log'):
                    file_path = os.path.join(library_analysis_path, file_name)
                    safe_remove(file_path)
            for subdir, dirs, files in os.walk(parent_dir):
                for dir in dirs:
                    subdir_path = os.path.join(subdir, dir)
                    analysis_dir = os.path.join(subdir_path, "FMOPhore")
                    if not os.path.isdir(analysis_dir):
                        return
                    for file in common_files_to_remove:
                        file_path = os.path.join(analysis_dir, file)
                        safe_remove(file_path)
                    pocketpqr = os.path.join(library_analysis_path, "pocket.pqr")
                    safe_remove(pocketpqr)
        parent_dir_path = "."            
        main(parent_dir_path)
if __name__ == "__main__":
    """
    FMOPhore V 0.1 - LibAnalysis - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    pdb_processor = LibAnalysis(self.pdb_file)
    pdb_processor.Lib_analyze()