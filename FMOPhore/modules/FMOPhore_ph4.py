import numpy as np
import sys
import os
import re
import glob
import shutil
import csv
import subprocess
import math
###############################################################################################################################
#                                           process binding site analysis
###############################################################################################################################
class Pharmacophores:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        
    def run_Pharmacophores(self):
        pdb_lines_ligand = self.load_processed_pdb_file()
        pdb_lines_ligand = self.rename_LIG(pdb_lines_ligand)
        self.save_processed_pdb_file(pdb_lines_ligand)
        self.process_BS_analysis()
        self.run_ph4()
    ################################################################
    def load_processed_pdb_file(self):
        pdb_file_path = f"{self.pdb_file[:-4]}.pdb"
        with open(pdb_file_path, 'r') as f:
            return f.readlines()
    ################################################################
    def rename_LIG(self, pdb_lines):
        excluded_strings = ['FMT', 'IOD', 'HOH', 'WAT', 'T3P', 'CL','BU3',
                            'MG', 'ZN','PEG', 'DMS', 'GOL', 'BO3', 
                            'EDO', 'SO4', 'ACE', 'NMA', 'NME', 'ACT', 
                            'MES', 'OCY','FMT','PEG','SEP', 'BME', 'CSO', 
                            'IMD','TPO', 'TCE', 'MDP', 'NI'
                            ]
        for i, line in enumerate(pdb_lines):
            if line.startswith('HETATM') and line[17:20].strip() not in excluded_strings:
                pdb_lines[i] = line[:22] + '9999' + line[26:]
        for i, line in enumerate(pdb_lines):
            record_type = line[:6].strip()
            if record_type in ['ATOM', 'HETATM']:
                if record_type == 'HETATM':
                    if line[17:20].strip() not in excluded_strings:
                        line = line[:17] + 'LIG' + line[20:]
                pdb_lines[i] = line
        return pdb_lines
    ################################################################
    def save_processed_pdb_file(self, pdb_lines):
        with open(f"{self.pdb_file[:-4]}.pdb", 'w') as f:
            f.writelines(pdb_lines)
    ################################################################
    #  process binding site analysis
    ################################################################
    def process_BS_analysis(self):
        cwd = os.getcwd()
        pdb_file = f"{self.pdb_file[:-4]}.pdb"
        list_files = subprocess.run(["fpocket", "-f", pdb_file])
        def trouve_ligands(input_filename, output_filename):
            ligands_set = set()
            with open(input_filename, "r") as pdb_file:
                for line in pdb_file:
                    if line.startswith('HETATM') and line[17:20].strip() != 'HOH':
                        ligands_set.add(f"{input_filename} {line[17:20].strip()}")
            unique_ligands = sorted(set(ligands_set))
            with open(output_filename, 'w') as output_file:
                output_file.write("\n".join(unique_ligands))
        trouve_ligands(pdb_file, self.pdb_file[:-4] + "_Ligands.txt")
        for filename in os.listdir(cwd):
            if filename.startswith(self.pdb_file[:-4] + "_Ligands.txt"):
                list_files = subprocess.run(["dpocket", "-f", filename, "-e", "-d 4.0"])
        def trouve_ligands_nom(fname):
            with open(fname, "r") as pdb_file:
                return [line.rstrip() for line in pdb_file if line[0:6].strip() == 'HETATM' and line[17:20].strip() != 'HOH']
        def process_pdb_files(directory):
            output_filename = self.pdb_file[:-4] + "_Ligands_APO.txt"
            with open(output_filename, 'w') as output_file:
                unique_lines_set = set()
                for root, dirs, files in os.walk(directory):
                    for ligand_file in [f for f in files if f.endswith("H.pdb") and f[-9:] != "lig_H.pdb"]:
                        protein_file = next((f for f in files if f.endswith("protein.pdb") and f[0:4] == ligand_file[0:4]), None)
                        if protein_file:
                            ligands_fp = trouve_ligands_nom(os.path.join(root, ligand_file))
                            for ligand_line in ligands_fp:
                                if ligand_line[0:6].strip() == 'HETATM' and ligand_line[17:20].strip() != "HOH":
                                    unique_line = f"{protein_file} {ligand_file} {ligand_line[17:20]}"
                                    unique_lines_set.add(unique_line)
                output_file.writelines(sorted(unique_lines_set))
        process_pdb_files(cwd)
        for root, subdirs, filenames in os.walk(cwd):
            for filename in filenames:
                if filename.startswith(self.pdb_file[:-4] + "_Ligands_APO.txt"):
                    list_files = subprocess.run(["tpocket", "-L", filename])
        os.remove(self.pdb_file[:-4] + "_Ligands_APO.txt")
        os.remove(self.pdb_file[:-4] + "_Ligands.txt")
        files_to_move = glob.glob('stats_p.txt') + glob.glob('Ph4_3D.txt') + glob.glob('stats_g.txt') + glob.glob('dpout_fpocketp.txt') + glob.glob('dpout_fpocketnp.txt') + glob.glob('dpout_explicitp.txt')
        for folder in glob.glob('*_out'):
            if folder.endswith('_out'):
                new_destination_folder = folder.replace(folder, 'FMOPhore')
                shutil.move(folder, new_destination_folder)
                break
        destination_folder = 'FMOPhore'
        for file_path in files_to_move:
            shutil.move(file_path, destination_folder)
        os.chdir('FMOPhore')
    # ################################################################################################################################
        file_paths = ['dpout_fpocketp.txt',
                     'dpout_fpocketnp.txt',
                     'dpout_explicitp.txt',
                     'stats_p.txt']
        for file_path in file_paths:
            with open(file_path, 'r') as file:
                lines = file.readlines()
            header = lines[0].strip().split()
            header = lines[0].strip().split()
            try:
                data = lines[1].strip().split()
                data_dict = {}
                for col_num, desc in enumerate(header[1:]):
                    data_dict[desc] = data[col_num + 1]
                max_desc_length = max(len(desc) for desc in data_dict)
                modified_descriptions = [f"{desc:{max_desc_length}}:   {data_dict[desc]}" for desc in data_dict]
                modified_content = "\n".join(modified_descriptions)
                modified_file_path = file_path.replace('.txt', '_modified.txt')
                with open(modified_file_path, 'w') as modified_file:
                    modified_file.write(modified_content)
                os.rename(modified_file_path, file_path)
            except IndexError:
                data = None
    # ################################################################################################################################
        def read_pdb(file_path):
            data = []
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        data.append(line)
            return data
        def find_matching_residues(pdb1_data, pdb2_data):
            matching_residues = []
            for line1 in pdb1_data:
                if line1[0:6] in ["ATOM  ", "HETATM"]:
                    res_id1 = line1[22:27].strip()
                    res_name1 = line1[17:20].strip()
                    for line2 in pdb2_data:
                        if line2[0:6] in ["ATOM  ", "HETATM"]:
                            res_id2 = line2[22:27].strip()
                            res_name2 = line2[17:20].strip()
                            if res_id1 == res_id2 and res_name1 == res_name2:
                                matching_residues.append((res_name1, res_id1))
            return matching_residues
        def is_druggable(score):
            return score > 0.3
        os.chdir("pockets")
        input_file_paths = glob.glob("*.pdb")
        pocket_scores = {}
        pocket_scores_file = f"../{self.pdb_file[:-4]}_info.txt"
        with open(pocket_scores_file, "r") as score_file:
            lines = score_file.readlines()
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith("Pocket"):
                    pocket_num = line.split()[1]
                    druggability_score = None
                    for j in range(i, len(lines)):
                        if "Score :" in lines[j]:
                            score_line = lines[j].strip().split(":")
                            if len(score_line) == 2:
                                druggability_score = float(score_line[1].strip())
                                break
                    if druggability_score is not None:
                        pocket_scores[pocket_num] = druggability_score
                i += 1
        concatenated_data = []
        for i, file_path in enumerate(input_file_paths):
            pdb1_data = read_pdb(file_path)
            pocket_num = file_path.split("_")[0][6:]  
            druggability_score = pocket_scores.get(pocket_num, 0.000)
            if is_druggable(druggability_score):
                matching_residues = []
                for j, other_file_path in enumerate(input_file_paths):
                    if i != j:
                        pdb2_data = read_pdb(other_file_path)
                        matching_residues.extend(find_matching_residues(pdb1_data, pdb2_data))
                concatenated_data.extend(pdb1_data)
                for line in pdb1_data:
                    res_id = line[22:27].strip()
                    res_name = line[17:20].strip()
                    if (res_name, res_id) in matching_residues:
                        concatenated_data.append(line)
        output_file_path = "../pocket.pdb"
        with open(output_file_path, "w") as file:
            seen_lines = set()  
            for line in concatenated_data:
                if line not in seen_lines:  
                    file.write(line)
                    seen_lines.add(line)
        input_file_paths = glob.glob("*.pqr")
        concatenated_data = []
        for i, file_path in enumerate(input_file_paths):
            pqr1_data = read_pdb(file_path)
            pocket_num = file_path.split("_")[0][6:]  
            druggability_score = pocket_scores.get(pocket_num, 0.000)
            if is_druggable(druggability_score):
                concatenated_data.extend(pqr1_data)
        output_file_path = "../pocket.pqr"
        with open(output_file_path, "w") as file:
            for line in concatenated_data:
                file.write(line)
        os.chdir("..")
        def distance(coord1, coord2):
            x1, y1, z1 = coord1
            x2, y2, z2 = coord2
            return ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**0.5
        def is_far(coord, all_coords, threshold=2):
            for other_coord in all_coords:
                if distance(coord, other_coord) <= threshold:
                    return False
            return True
        all_growth_vectors = 'pocket.pqr'
        pocket = 'pocket.pdb'
        ligand = f"../{self.pdb_file[:-4]}_FMO_lig_H.pdb"
        single_growth_vectors = 'single_growth_vectors.pdb'
        single_growth_vectors_volume = 'single_growth_vectors_volume.txt'
        pocket_coords = []
        with open(pocket, 'r') as file1:
            for line in file1:
                if line.startswith('ATOM'):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    pocket_coords.append((x, y, z))
        file1_coords = []
        with open(all_growth_vectors, 'r') as file1:
            for line in file1:
                if line.startswith('ATOM'):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    file1_coords.append((x, y, z))
        file2_coords = []
        with open(ligand, 'r') as file2:
            for line in file2:
                if line.startswith('HETATM'):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    file2_coords.append((x, y, z))
        filtered_coords = []
        renumber = 1 
        for line in open(all_growth_vectors, 'r'):
            if line.startswith('ATOM'):
                line = line[:22] + str(renumber).rjust(4) + line[26:]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coord = (x, y, z)
                if is_far(coord, file2_coords):
                    filtered_coords.append(line)
                    renumber += 1  
        with open(single_growth_vectors, 'w') as output_file:
            output_file.writelines(filtered_coords)
        def center_of_mass(coords):
            num_coords = len(coords)
            if num_coords == 0:
                return None
            sum_x = sum(y for x, y, z in coords)
            sum_y = sum(y for x, y, z in coords)
            sum_z = sum(z for x, y, z in coords)
            center_x = sum_x / num_coords
            center_y = sum_y / num_coords
            center_z = sum_z / num_coords
            return center_x, center_y, center_z
        def sphere_volume(radius):
            return (4/3) * math.pi * (radius ** 3)
        total_volume = 0
        for line in filtered_coords:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords = [(x, y, z)]  
            center = center_of_mass(coords)
            if center is not None:
                sphere_vol = sphere_volume(1.2)
                total_volume += sphere_vol
        with open(single_growth_vectors_volume, 'w') as volume_file:
            volume_file.write(str(total_volume))
        with open(pocket, "r") as f1:
            lines1 = f1.readlines()
        lines1 = [line for line in lines1 if 'HEADER' not in line]
        lines1 = [line for line in lines1 if 'END' not in line]
        lines1 = [line for line in lines1 if 'TER' not in line]
        with open(single_growth_vectors, "r") as f3:
            lines3 = f3.readlines()
        lines3 = [line for line in lines3 if 'HEADER' not in line]
        lines3 = [line for line in lines3 if 'END' not in line]
        lines3 = [line for line in lines3 if 'TER' not in line]
        def extract_coords(line):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            return np.array([x, y, z])
        coords1 = [extract_coords(line) for line in lines1 if line.startswith("ATOM")]
        coords3 = [extract_coords(line) for line in lines3 if line.startswith("ATOM")]
        def calculate_center_of_mass(coords):
            total_atoms = len(coords)
            center_of_mass = np.sum(coords, axis=0) / total_atoms
            return center_of_mass
        center_of_mass1 = calculate_center_of_mass(coords1)
        for atom_id, coords3 in enumerate(coords3, start=1):
            center_of_mass3 = calculate_center_of_mass([coords3])
            nearest_id = min(range(len(coords1)), key=lambda k: np.linalg.norm(coords1[k] - center_of_mass3))
            new_line = lines3[atom_id - 1][:6] + lines1[nearest_id][6:26] + lines3[atom_id - 1][26:]
            lines3[atom_id - 1] = new_line
        with open(single_growth_vectors, "w") as f3_updated:
            f3_updated.writelines(lines3)
        with open(single_growth_vectors, "r") as f3_updated:
            lines_updated = f3_updated.readlines()
        modified_lines = []
        for line in lines_updated:
            record_type = line[0:6].strip()
            atom_id = line[6:11].strip()
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain_id = line[21]
            residue_number = line[22:26].strip()
            insertion_code = line[26]
            x_coord = line[30:38]
            y_coord = line[38:46]
            z_coord = line[46:54]
            occupancy = line[54:60]
            temp_factor = line[60:66]
            element = line[76:78].strip()
            charge = line[78:80].strip()
            if record_type == "ATOM" or record_type == "HETATM":
                if atom_name.startswith("N"):
                    atom_name = "O"
                elif atom_name.startswith("O"):
                    atom_name = "N"
                elif atom_name.startswith("C"):
                    if residue_name not in ["HIS", "TYR", "TRP", "PHE"]:
                        atom_name = "C"
                if residue_name in ["HIS", "TYR", "TRP", "PHE"] and atom_name.startswith("C"):
                    atom_name = "Pi"
            formatted_line = f"{record_type:<6}{atom_id:>5}{atom_name:>5}{residue_name:>4} {chain_id}{residue_number:>4}{insertion_code}   {x_coord}{y_coord}{z_coord}{occupancy}{temp_factor}          {element:<2}{charge:2}\n"
            modified_lines.append(formatted_line)
        sorted_lines = sorted(modified_lines, key=lambda line: int(line[22:26]))
        with open(single_growth_vectors, "w") as f3_modified:
            f3_modified.writelines(sorted_lines)
        os.remove(f"{self.pdb_file[:-4]}" + "_pockets.pqr")
        os.remove(f"{self.pdb_file[:-4]}" + "_PYMOL.sh")
        os.remove(f"{self.pdb_file[:-4]}" + "_out.pdb")
        os.remove(f"{self.pdb_file[:-4]}" + "_VMD.sh")
        os.remove(f"{self.pdb_file[:-4]}" + ".pml")
        os.remove(f"{self.pdb_file[:-4]}" + ".tcl")
        for file_name in glob.glob('pockets/*.pqr'):
            os.remove(file_name)
    ################################################################
    def run_ph4(self):
        os.system(f'plip -f ../{self.pdb_file[:-4]}.pdb -t -q -s > .log')
        os.remove(f'{self.pdb_file[:-4]}_protonated.pdb')
        os.remove(f'.log')
        for file_name in os.listdir('.'):
            if file_name.startswith('plip') and file_name.endswith('.pdb'):
                file_path = os.path.join('.', file_name)
                os.remove(file_path)
        with open("report.txt", "r") as file:
            rows = []
            keyword = None  
            section_title = None  
            section_title_end = None  
            for line in file:
                if re.match(r'^[+=]+\n$', line) or re.match(r'^[+-]+\n$', line):
                    continue
                csv_line = line.replace('|', '|').replace('LIGCOO', 'LIGCOO.x,LIGCOO.y,LIGCOO.z').replace('PROTCOO', 'PROTCOO.x,PROTCOO.y,PROTCOO.z')
                csv_line = re.sub(r'\s*,\s*', ',', csv_line)
                csv_line = csv_line.strip(',')
                if csv_line.startswith('RESNR'):
                    csv_line += ',Bond_Type'
                rows.append(csv_line)
                if not section_title and line.startswith('**'):
                    section_title = re.search(r'\*\*(.*?)\*\*', line)
                    if section_title:
                        section_title = section_title.group(1)
            output_lines = []
            for i, row in enumerate(rows):
                output_line = row.strip()
                if section_title and not row.startswith('**') and section_title_end is not None and i > section_title_end:
                    output_line += ',' + section_title
                output_lines.append(output_line)
                if i == section_title_end:
                    output_lines.append('')
        output_lines = [line if not line.startswith('""') else '' for line in output_lines]
        with open("Report.csv", "w", newline='') as outfile:
            writer = csv.writer(outfile)
            for row in output_lines:
                writer.writerow(row.split(','))
        os.remove(f'report.txt')
        with open("Report.csv", "r") as infile:
            lines = infile.readlines()
        output_lines = []
        for line in lines:
            if not line.startswith('""') and "Prediction of noncovalent interactions for PDB structure" not in line and "Created on" not in line and "If you are using PLIP in your work,please cite:" not in line and "Analysis was done on model" not in line and "Adasme,M. et al. PLIP 2021" not in line and "Interacting chain(s)" not in line:
                output_lines.append(line)
            else:
                output_lines.append('\n')
        with open("Report.csv", "w") as outfile:
            outfile.writelines(output_lines)
        with open("Report.csv", "r") as infile:
            lines = infile.readlines()
        output_lines = []
        add_keyword = False
        keyword = ""
        for line in lines:
            if line.startswith("**"):
                keyword = line.strip()[2:-2]
                keyword = "_".join(keyword.split())
                add_keyword = True
                output_lines.append(line)
            elif line.strip() != "":
                if add_keyword:
                    line = line.strip() + keyword + "\n"
                output_lines.append(line)
            else:
                output_lines.append('\n')
        non_empty_lines = [line.strip() for line in output_lines if line.strip()]
        exclusion_keywords = {
            'Excluded molecules as ligands:', 'No interactions detected.', 'DMS', 'BO3', 'GOL',
            'AGOL', 'BGOL', 'ZN', 'EDO', 'SO4', 'ACE', 'NMA', 'NME', 'ACT', 'MES', 'OCY',
            'SMALLMOLECULE'}
        non_empty_lines = [
            line for line in non_empty_lines
            if not any(line.startswith(kw) for kw in ['**', 'Excluded molecules as ligands:', 'No interactions detected.'])
            and 'LIG' in line
            and line[16:17] not in {'B', 'C', 'D', 'E'}
            and not any(kw in line for kw in exclusion_keywords)]
        modified_lines = []
        for line in non_empty_lines:
            if line.startswith("RESNR"):
                parts = line.split("|")
                parts[-1] = ""
                parts[4] = "LIG"
                parts[3] = parts[4]
                parts[4] = ""
                parts[0] = parts[0] + "-" + parts[1]
                parts[1] = ""
                parts[2] = parts[2]
                modified_lines.append("|".join(parts).rstrip("|").replace("||", "|"))
            elif line[2].isdigit():
                parts = line.split("|")
                parts[0] = parts[2] + parts[3] + parts[1]
                modified_lines.append("|".join(parts).replace("||", "|"))
            else:
                modified_lines.append(line)
        with open("Report.csv", "w") as outfile:
            outfile.write('\n'.join(modified_lines))
        with open("Report.csv", "r") as infile:
            lines = infile.readlines()
        modified_lines = []
        lines = [line for line in lines if 'RESNR' not in line]
        for line in lines:
            if not line.startswith("RESNR"):
                parts = line.split("|")
                modified_line = ""
                if 'Hydrophobic_Interactions' in line:
                    modified_line = parts[0]+"\t"+parts[7]+"\t"+"0"+"\t"+"\t".join(parts[-3:])
                elif 'Hydrogen_Bonds' in line:
                    modified_line = parts[0]+"\t"+parts[8]+"\t"+parts[10]+"\t"+"\t".join(parts[-3:])
                elif 'pi-Stacking' in line:
                    modified_line = parts[0]+"\t"+parts[8]+"\t"+parts[9]+"\t"+"\t".join(parts[-3:])
                elif 'pi-Cation_Interactions' in line:
                    modified_line = parts[0]+"\t"+parts[8]+"\t"+"0"+"\t"+"\t".join(parts[-3:])
                elif 'Water_Bridges' in line:
                    modified_line = parts[0]+"\t"+parts[7]+"\t"+parts[9]+"\t"+parts[-4]+"\t"+parts[-3]+"\t"+parts[-1]
                elif 'Salt_Bridges' in line:
                    modified_line = parts[0]+"\t"+parts[8]+"\t"+"0"+"\t"+"\t".join(parts[-3:])
                elif 'Halogen_Bonds' in line:
                    modified_line = parts[0]+"\t"+parts[8]+"\t"+parts[9]+"\t"+parts[-2]+"\t"+parts[-3]+"\t"+parts[-1]
                elif 'Metal_Complexes' in line:
                    modified_line = parts[0] + "\t" + parts[9]  + "\t" + "0"       + "\t"  + "\t".join(parts[-7:])
                modified_lines.append(modified_line.strip("\n"))
            else:
                modified_lines.append(line.strip("\n"))
        os.remove(f'Report.csv')
        with open("data.txt", "w") as outfile:
            outfile.write("Fragment\tdistance\tAngle\tLIG.x\tLIG.y\tLIG.z\tPROT.x\tPROT.y\tPROT.z\tBond_type\n")
            modified = []
            for line in modified_lines:
                parts = line.split()
                modified_line = "".join(parts[:3]) + '\t' + '\t'.join(parts[3:])
                modified_line = modified_line.replace(' ', '\t').replace(',', '\t')
                modified.append(modified_line)
            outfile.write('\n'.join(modified))
        ###################################################################################################
        def extract_data_from_pocket(file_path):
            extracted_data_protein = []
            extracted_data_ligand = []
            extracted_data_ligand_atom_index = []
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith("ATOM"):
                        space = " " in line[line.index("ATOM") + 22 : line.index("ATOM") + 23]
                        if space:
                            match = re.search(r"ATOM\s+\d+\s+(\S+)\s+(\w+)\s+(\w+)\s+(\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)", line)
                            if match:
                                residue = match.group(2)
                                atom_name = match.group(1)
                                chain_id = match.group(3)
                                residue_id = match.group(4)
                                x_coord = match.group(5)
                                y_coord = match.group(6)
                                z_coord = match.group(7)
                                extracted_data_protein.append(f"{residue}{chain_id}{residue_id} {atom_name} {float(x_coord):8.3f}{float(y_coord):8.3f}{float(z_coord):8.3f}")
                        elif not space:
                            match = re.search(r"ATOM\s+\d+\s+(\S+)\s+(\w+)\s+((\w+)\d+)+\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)", line)
                            if match:
                                residue = match.group(2)
                                atom_name = match.group(1)
                                chain_id = match.group(3)
                                residue_id = match.group(4)
                                x_coord = match.group(5)
                                y_coord = match.group(6)
                                z_coord = match.group(7)
                                extracted_data_protein.append(f"{residue}{chain_id}{residue_id} {atom_name} {float(x_coord):8.3f}{float(y_coord):8.3f}{float(z_coord):8.3f}")
                    elif line.startswith("HETATM"):
                        space = " " in line[line.index("HETATM") + 6 : line.index("HETATM") + 7]
                        if space:
                            match = re.search(r"HETATM\s+\d+\s+(\S+)\s+(\w+)\s+((\w+)\d+)+\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)", line)
                            if match and all(match.group(i) is not None for i in range(1, 7)):
                                residue = match.group(2)
                                atom_name = match.group(1)
                                chain_id = match.group(3)
                                residue_id = match.group(4)
                                x_coord = match.group(5)
                                y_coord = match.group(6)
                                z_coord = match.group(7)
                                extracted_data_ligand.append(f"{residue}{chain_id}{residue_id} {atom_name} {float(x_coord):8.3f}{float(y_coord):8.3f}{float(z_coord):8.3f}")
                        elif not space:
                            match = re.search(r"HETATM\d+\s+(\S+)\s+(\w+)\s+((\w+)\d+)+\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)", line)
                            if match and all(match.group(i) is not None for i in range(1, 7)):
                                residue = match.group(2)
                                atom_name = match.group(1)
                                chain_id = match.group(3)
                                residue_id = match.group(4)
                                x_coord = match.group(5)
                                y_coord = match.group(6)
                                z_coord = match.group(7)
                                extracted_data_ligand.append(f"{residue}{chain_id}{residue_id} {atom_name} {float(x_coord):8.3f}{float(y_coord):8.3f}{float(z_coord):8.3f}")
            return extracted_data_protein + extracted_data_ligand
        ###################################################################################################
        def compare_and_update_protdata(data, data_file_path):
            updated_data_prot = []
            with open(data_file_path, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    columns = line.split()
                    residue = columns[0]
                    PROTx = columns[-4]
                    PROTy = columns[-3]
                    PROTz = columns[-2]
                    for data_line in data:
                        data_columns = data_line.split()
                        if data_columns[0].startswith(residue) and data_columns[-3] == PROTx and data_columns[-2] == PROTy and data_columns[-1] == PROTz:
                            updated_line = line.strip() + f"\t{data_columns[1]}\n"
                            updated_data_prot.append(updated_line)
                            break
                    else:
                        updated_data_prot.append(line)
            return updated_data_prot
        def compare_and_update_ligdata(data, data_file_path):
            updated_data_lig = []
            with open(data_file_path, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    columns = line.split()
                    residue = columns[0]
                    LIGx = columns[-8]
                    LIGy = columns[-7]
                    LIGz = columns[-6]
                    for data_line in data:
                        data_columns = data_line.split()
                        if data_columns[-3] == LIGx and data_columns[-2] == LIGy and data_columns[-1] == LIGz:
                            updated_line = line.strip() + f"\t{data_columns[1]}\n"
                            updated_data_lig.append(updated_line)
                            break
                    else:
                        updated_data_lig.append(line)
            return updated_data_lig
        file_path = f"../{self.pdb_file[:-4]}.pdb" 
        extracted_data_combined = extract_data_from_pocket(file_path)
        data_file_path = "data.txt"  
        updated_data_prot = compare_and_update_protdata(extracted_data_combined, data_file_path)
        updated_data_prot[0] = updated_data_prot[0].strip() + "\tprot_atom_id\n"
        with open(data_file_path, 'w') as file:
            file.writelines(updated_data_prot)
        input_file_path = data_file_path  
        output_file_path = "modified_file.txt"  
        with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
            for line in input_file:
                if "pi-Stacking" in line:
                    line = line.strip() + "\tPi\n"
                elif "Salt_Bridges" in line and "N" not in line:
                    line = line.strip() + "\tSa\n"
                elif "Salt_Bridges" in line and "N" in line:
                    line = line.strip() + "\tO\n"
                elif "pi-Cation_Interactions" in line and "N" not in line:
                    line = line.strip() + "\tPi\n"
                output_file.write(line)
        data_file_path = "modified_file.txt"  
        updated_data_lig = compare_and_update_ligdata(extracted_data_combined, output_file_path)
        updated_data_lig[0] = updated_data_lig[0].strip() + "\tlig_atom_id\n"
        with open(data_file_path, 'w') as file:
            file.writelines(updated_data_lig)
        input_file_path = data_file_path  
        output_file_path = "Ph4_3D.txt"  
        c_with_number_pattern = re.compile(r'C\d+')
        with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
            for line in input_file:
                if "pi-Stacking" in line:
                    line = line.strip() + "\tPi\n"
                elif "pi-Cation_Interactions" in line and ("Pi" not in line):
                    line = line.strip() + "\tPi\n"
                elif "pi-Cation_Interactions" in line and ("Pi" in line and "N" not in line and not c_with_number_pattern.search(line)):
                    line = line.strip() + "\tN\n"
                elif "Salt_Bridges" in line and line.strip().endswith("Sa"):
                    line = line.strip() + "\tO\n"
                output_file.write(line)
        os.remove(f"data.txt")
        os.remove(f"modified_file.txt")
        ###################################################################################################
        def process_file_to_pdb(input_file_path, output_file_path):
            atom_mapping = {
                'Pi': 'R', 
                'O': 'A',
                'C': 'H',
                'N': 'D',
                'NH': 'P',
                'O': 'N'
            }
            with open(input_file_path, 'r') as input_file:
                lines = input_file.readlines()
            with open(output_file_path, 'w') as output_file:
                output_file.write(f"{output_file_path}\n")
                for line in lines[1:]:
                    columns = line.strip().split('\t')
                    original_name = ''.join(char for char in columns[-1].strip()[0:2] if char.isalpha())
                    if original_name in atom_mapping:
                        name = atom_mapping[original_name]
                    else:
                        name = original_name
                    coordinates = columns[6:9]
                    output_file.write(f"{name}, {','.join(coordinates)}\n")
        input_file_path = 'Ph4_3D.txt'  
        output_file_path = f'{self.pdb_file[:-4]}_ph4_3D.txt'  
        process_file_to_pdb(input_file_path, output_file_path)
        custom_order = {'R': 0, 'D': 1, 'A': 2, 'P': 3, 'N': 4, 'H': 5}
        with open(output_file_path, 'r') as file:
            lines = [line.strip() for line in file if line.strip()]
        def is_valid(line):
            parts = line.split(',')
            return len(parts) == 4 and parts[0].strip() in custom_order
        sorted_lines = sorted([line for line in lines if is_valid(line)], key=lambda line: (custom_order[line.split(',')[0]], [float(coord.strip()) for coord in line.split(',')[1:]]))
        with open(output_file_path, 'w') as file:
            file.write(f"{output_file_path}\n")
            file.write('\n'.join(sorted_lines))
################################################################################################################################
if __name__ == "__main__":
    """
    FMOPhore V.0.1 - Pharmacophores - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    pdb_processor = Pharmacophores(self.pdb_file)
    pdb_processor.run_Pharmacophores()