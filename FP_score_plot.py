import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import numpy as np
import matplotlib.cm as cm
import os

def format_fragment_name(fragment):
    residue = fragment[:3].capitalize()
    chain = fragment[3]
    number = fragment[4:]
    return f"{residue}-{chain}-{number}"

data = pd.read_csv('merged_file.csv')
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
        texts.append(ax1.text(row['frequency'], row['Energy'], row['Formatted_Fragment'], fontsize=14))  # Fontsize reduced for text

adjust_text(texts, ax=ax1, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))
ax1.axhline(mean_energy, color='grey', linestyle='--')
ax1.axvline(mean_frequency, color='grey', linestyle='--')
ax1.set_xlabel('Normalization of Percentage of interaction (PI%)', fontsize=14, fontweight='bold', color='black')  # Reduced fontsize
ax1.set_ylabel('Normalization of Interaction Energy (ΔE)', fontsize=14, fontweight='bold', color='black')  # Reduced fontsize

ax1.tick_params(axis='x', labelsize=12)  
ax1.tick_params(axis='y', labelsize=12)
plt.setp(ax1.get_xticklabels(), fontweight='bold')
plt.setp(ax1.get_yticklabels(), fontweight='bold')
plt.tight_layout()
plt.savefig('FP_score_vs_Hotspots.png', dpi=300)  


def format_fragment_name(fragment):
    residue = fragment[:3].capitalize()
    chain = fragment[3]
    number = fragment[4:]
    return f"{residue}-{chain}-{number}"

file_path = 'merged_file.csv'  
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

for bond_type in bond_types:
    if bond_type not in df.columns:
        df[bond_type] = 0

for i, row in df.iterrows():
    total_bond_interactions = row[bond_types].sum()
    if total_bond_interactions > 0:
        df.loc[i, bond_types] = (row[bond_types] / total_bond_interactions) * row["Percentage of Interaction (%)"]

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

ax1.set_ylim(0, 100)
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
output_file_coolwarm = 'Binding_site_residues_barchart_normalized.png'
plt.savefig(output_file_coolwarm, dpi=400) 

def format_fragment_name(fragment):
    residue = fragment[:3].capitalize()  
    chain = fragment[3]  
    number = fragment[4:] 
    return f"{residue}-{chain}-{number}"

parent_dir = os.getcwd()
selected_data_path = os.path.join(parent_dir, 'selected_data.csv')
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
plt.figure(figsize=(12, 6))  # Reduced figure size

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
heatmap_filename = os.path.join(parent_dir, 'Ph4_heatmap.png')
plt.tight_layout()
plt.savefig(heatmap_filename, dpi=400)  

def get_residue_number(formatted_fragment):
    return int(formatted_fragment.split('-')[-1])

data = pd.read_csv('merged_file.csv')
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
df['Residue_Number'] = df['Formatted_Fragment'].apply(get_residue_number)
df_sorted = df.sort_values(by='Residue_Number')
unique_fragments_sorted = df_sorted['Formatted_Fragment'].unique()
palette = sns.color_palette("hsv", len(unique_fragments_sorted))
color_map = dict(zip(unique_fragments_sorted, palette))
fig2, ax2 = plt.subplots(figsize=(8, 5))
colors = [color_map[fragment] for fragment in df_sorted['Formatted_Fragment']]
ax2.bar(df_sorted['Formatted_Fragment'], df_sorted['FP_score'], width=0.4, alpha=0.6, color=colors, edgecolor='black')
ax2.axhline(threshold, color='black', linestyle='--', label=f'Threshold: {threshold:.2f}')
ax2.set_xlabel('Binding site Hotspots', fontsize=14, fontweight='bold')
ax2.set_ylabel('FP_score', fontsize=14, fontweight='bold', color='black')
ax2.tick_params(axis='y', labelcolor='black', labelsize=12)
ax2.tick_params(axis='x', labelsize=12, rotation=90)
ax2.set_title('FP-score vs Binding site Hotspots', fontsize=16, fontweight='bold')
plt.setp(ax2.get_xticklabels(), fontweight='bold')
plt.setp(ax2.get_yticklabels(), fontweight='bold')
plt.tight_layout()
plt.savefig('Energy_vs_Percentage_of_Interaction.png', dpi=300)
print("Plot saved successfully.")

fp_data = pd.read_csv('FP_scores.csv')
mean_frequency = fp_data['frequency'].mean()
mean_energy = fp_data['Energy'].mean()
filtered_fp_data = fp_data[(fp_data['frequency'] > mean_frequency) | (fp_data['Energy'] > mean_energy)]
fragments_to_keep = set(filtered_fp_data['Fragment'].apply(format_fragment_name))
print(fragments_to_keep)
pdb_file = 'Ph4_lig.pdb'  
output_pdb_file = 'hotspot_cloud.pdb'
with open(pdb_file, 'r') as infile, open(output_pdb_file, 'w') as outfile:
    for line in infile:
        if line.startswith('ATOM'):
            residue = line[17:20].strip()
            chain = line[21].strip()
            number = line[22:26].strip()
            fragment_name = format_fragment_name(f"{residue}{chain}{number}")
            if fragment_name in fragments_to_keep:
                outfile.write(line)

print(f"Filtered PDB file saved to: {output_pdb_file}")
