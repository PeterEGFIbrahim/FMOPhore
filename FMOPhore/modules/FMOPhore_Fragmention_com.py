import sys
import os
import argparse
import subprocess
import shutil, os, glob

###################################################################################################
class Fragmention_Processor:
    def __init__(self, input_pdb_file, output_file='output.inp', qm=None):
        self.input_pdb_file = input_pdb_file
        self.output_file = output_file
        self.input_pdb_lines = [] 
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.qm = qm
    def process_Fragmention(self):
        qm_file_path = os.path.join(self.script_dir, f'{self.qm}.md')
        with open(qm_file_path, 'r') as qm_file:
            content = qm_file.read()
        sections = content.split('\n\n')
        header_content, body_content, tail_content = sections[0], sections[1], sections[2]
        chain_data = self.read_pdb_file()
        for chain_id, chain_lines in chain_data.items():
            chain_dir = os.path.join('chains', chain_id)
            os.makedirs(chain_dir, exist_ok=True)
            self.process_pdb_file(chain_lines, chain_dir)
        output_directory = 'chains'
        self.concatenate_cord_files_sorted(output_directory)
        for chain_id in chain_data.keys():
            chain_dir = os.path.join('chains', chain_id)
            shutil.rmtree(chain_dir)
        self.generate_output_file(header_content, body_content, tail_content)
    def read_pdb_file(self):
        chains = {}
        current_chain = None
        with open(self.input_pdb_file, 'r') as f:
            self.input_pdb_lines = f.readlines()  
            for line in self.input_pdb_lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain_id = line[21]
                    if chain_id != current_chain:
                        current_chain = chain_id
                        if current_chain not in chains:
                            chains[current_chain] = []
                    chains[current_chain].append(line)
        return chains
    def process_pdb_file(self, chain_lines, chain_dir):
        os.makedirs(chain_dir, exist_ok=True)
        output_file_name = os.path.join(chain_dir, 'output.pdb')
        with open(output_file_name, 'w') as f:
            for i, line in enumerate(self.input_pdb_lines):
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    res_name = line[17:20].strip()
                    res_num = line[22:26].strip()
                    if res_name == 'ACE':
                        for j in range(i+1, len(self.input_pdb_lines)):
                            next_line = self.input_pdb_lines[j]
                            next_res_name = next_line[17:20].strip()
                            if next_res_name != 'ACE':
                                next_res_num = next_line[22:26].strip()
                                new_line = line[:17] + next_res_name + line[20:23] + next_res_num + line[26:]
                                f.write('{:<17s}{:>3s}{:1s}{:>3s}{}'.format(line[:17], next_res_name, line[20:23], next_res_num, line[26:-1]) + '\n')
                                break
                    else:
                        f.write(line)
                else:
                    f.write(line)
        with open(output_file_name, 'r+') as f:
            pdb_lines = f.readlines()
            f.seek(0)
            for i, line in enumerate(pdb_lines):
                if line.startswith('HETATM') and line[17:20] != 'LIG':
                    line = '{:<6s}{}'.format('ATOM', line[6:])
                    pdb_lines[i] = line
                f.write(line)
            f.truncate()
            f.seek(0)
            f.writelines(pdb_lines)
        pdb_file = open(output_file_name, 'r+')
        pdb_lines = pdb_file.readlines()
        prev_res_name = None
        prev_res_num = None
        pdb_file.seek(0)
        for line in pdb_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                res_name = line[17:20].strip()
                res_num = line[22:26].strip()
                if res_name == 'NMA':
                    if prev_res_name is not None and prev_res_num is not None:
                        new_line = line.replace('NMA', prev_res_name).replace(res_num, prev_res_num)
                        new_line = new_line[:26] + new_line[26:].replace('A', ' ')
                        pdb_file.write(new_line)
                        pdb_file.write(new_line)
                    else:
                        pdb_file.write(line)
                else:
                    pdb_file.write(line)
                    prev_res_name = res_name
                    prev_res_num = res_num
            else:
                if line.startswith('TER'):
                    line = line[:3] + '\n'
                pdb_file.write(line)
        pdb_file.truncate()
        pdb_file.close()
        unorganized_file_path = os.path.join(chain_dir, "unorganized.txt")
        with open(unorganized_file_path, "w") as f:
            sys.stdout = f
            with open(output_file_name, "r") as f:
                lines = f.readlines()
            first_atom_dict = {}
            last_atom_dict = {}
            for i, line in enumerate(lines):
                if line.startswith("ATOM"):
                    if line[12:16].strip() in ["N", "CA", "C", "O"]:
                        if (i+3<len(lines) and
                            lines[i+1][12:16].strip() == "CA" and
                            lines[i+2][12:16].strip() == "C" and
                            lines[i+3][12:16].strip() == "O"):
                            n_number = int(lines[i][6:11].strip())
                            ca_number = int(lines[i+1][6:11].strip())
                            c_number = int(lines[i+2][6:11].strip())
                            o_number = int(lines[i+3][6:11].strip())
                            res_name = line[17:20].strip()
                            res_num = line[22:26].strip()
                            res_key = f"{res_name}:{res_num}"
                            if res_key not in first_atom_dict:
                                first_atom_dict[res_key] = n_number
                            last_atom_dict[res_key] = o_number
            for key, value in first_atom_dict.items():
                print(f"{key}:{value}:N")
                print(f"{key}:{value+1}:CA")
                print(f"{key}:{value+2}:C")
                print(f"{key}:{value+3}:O")
                print(f"{key}:{value+4}:after")
            with open(output_file_name, "r") as f:
                lines = f.readlines()
            res_num = ""
            res_name = ""
            first_atom_dict = {}
            last_atom_dict = {}
            backbone_atom_dict = {}
            for line in lines:
                if line.startswith("ATOM"):
                    if (line[17:20].strip() == res_name and line[22:26].strip() == res_num and
                        line[11:16].strip() in ["N", "CA", "C", "O"]):
                        atom_num = int(line[6:11].strip())
                        backbone_atom_dict[line[11:16].strip()] = atom_num
                    res_name = line[17:20].strip()
                    res_num = line[22:26].strip()
                    curr_res_key = f"{res_name}:{res_num}"
                    if curr_res_key not in first_atom_dict:
                        first_atom_dict[curr_res_key] = int(line[6:11].strip())
                    last_atom_dict[curr_res_key] = int(line[6:11].strip())
            for key, value in first_atom_dict.items():
                print(f"{key}:{value}:first")
                print(f"{key}:{last_atom_dict[key]}:last")
                print()
            sys.stdout = sys.__stdout__
        ###################################################################################################
        with open(unorganized_file_path, 'r') as f:
            input_list = f.readlines()
        sorted_list = sorted(input_list, key=lambda x: int(x.split(":")[2]) if x.count(":") >= 3 else float('inf'))
        organized_file_path = os.path.join(chain_dir, "organized.txt")
        with open(organized_file_path, 'w') as f:
            for line in sorted_list:
                f.write(line)
        ###################################################################################################
        sorted_file_path = os.path.join(chain_dir, "sorted.txt")
        with open(sorted_file_path, "w") as f:
            sys.stdout = f
            flag_order = ["first", "N", "CA", "C", "O", "after", "last"]
            residue_dict = {}
            with open(organized_file_path, "r") as f:
                lines = f.readlines()
            lines = [line for line in lines if line.strip()]
            for line in lines:
                fields = line.strip().split(":")
                residue_id = fields[0]
                try:
                    residue_index = int(fields[1])
                except ValueError:
                    try:
                        residue_index = int(fields[1].lstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
                    except ValueError:
                        continue
                flag = fields[3]
                value = fields[2]
                if residue_index not in residue_dict:
                    residue_dict[residue_index] = {}
                residue_dict[residue_index][flag] = value
            for residue_index in sorted(residue_dict.keys()):
                residue_info = residue_dict[residue_index]
                output1 = []
                output2 = []
                for flag in flag_order:
                    if flag in residue_info:
                        if flag == "C":
                            output2.append(residue_info[flag])
                        elif flag == "O":
                            output2.append(residue_info[flag])
                        else:
                            output1.append(residue_info[flag])
                print(" ".join(output1))
                print(" ".join(output2))
            sys.stdout = sys.__stdout__
        ###################################################################################################
        coords_file_path = os.path.join(chain_dir, "coords.txt")
        with open(coords_file_path, "w") as f:
            sys.stdout = f
            filename = sorted_file_path
            with open(sorted_file_path, 'r') as f:
                lines = f.readlines()
            output = ""
            skip_next = False 
            for i in range(len(lines)):
                if skip_next:
                    skip_next = False
                    continue
                if len(lines[i].split()) == 2 and i+1 < len(lines):
                    output += lines[i].strip() + " " + lines[i+1]
                    skip_next = True
                else:
                    output += lines[i]
            sys.stdout = sys.__stdout__
        ###################################################################################################
        data = {}
        with open(organized_file_path, "r") as f:
            lines = [line for line in f if line.strip()]
            for line in lines:
                line = line.strip()
                parts = line.split(":")
                key = parts[0] + ":" + parts[1]
                value = parts[2]
                if key in data:
                    data[key].append(value)
                else:
                    data[key] = [value]
        modified_lines = [] 
        for key in data:
            nums = data[key] 
            if any(key.startswith(prefix) for prefix in ('HOH:', 'ZN:', 'MG:', 'CL:', 'NI:', 'CO:')):
                modified_line = f"{key}:{':'.join(nums)}:0:different" 
            else:
                flag = "same" if nums[0] == nums[1] else "different"  
                nums.pop(1)  # remove the second number
                modified_line = key + ":" + ":".join(nums) + ":0" + ":" + flag   
            modified_lines.append(modified_line) 
        modified_files0 = os.path.join(chain_dir,"modified_file0.txt")
        with open(modified_files0, "w") as new_file:
            new_file.write("\n".join(modified_lines))
        ###################################################################################################
        with open(modified_files0, "r") as f:
            lines = f.readlines()
        for i in range(len(lines)-1):
            line = lines[i].strip().split(":")
            residue_name = line[0]
            try:
                residue_number = int(line[1])
            except ValueError:
                try:
                    residue_number = int(line[1].lstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
                except ValueError:
                    continue
            if i == 0:
                new_residue_number = residue_number 
                line[-1] = "different"
                line[1] = str(new_residue_number)
                lines[i] = ":".join(line) + "\n"
            else:
                try:
                    next_residue_number = int(lines[i+1].strip().split(":")[1])
                except ValueError:
                    try:
                        next_residue_number = int(lines[i+1].strip().split(":")[1].lstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
                    except ValueError:
                        continue 
                if residue_number + 1 != next_residue_number:
                    lines[i+1] = lines[i+1].strip().split(":")
                    lines[i+1][-1] = "different"
                    lines[i+1] = ":".join(lines[i+1]) + "\n"
        modified_files1 = os.path.join(chain_dir,"modified_file1.txt")
        with open(modified_files1, "w") as f:
            f.writelines(lines)
        ###################################################################################################
        with open(modified_files1, "r") as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            parts = line.strip().split(":")
            residue_info = parts[2:]
            new_line = ",".join(residue_info)
            new_lines.append(new_line)
        modified_files2 = os.path.join(chain_dir, "modified_file2.txt")
        with open(modified_files2, "w") as f:
            f.write("\n".join(new_lines))
        ##################################################################################################
        with open(modified_files2, 'r') as f:
            lines = f.readlines()
        output_lines = []
        for i in range(len(lines)-1):
            line = lines[i]
            next_line = lines[i+1]
            if ',different' in line and ',same' in next_line:
                nums = line.split(',')
                third_num = nums[2]
                fourth_num = nums[3]
                output_lines.append(line)
                output_lines.append(f"{third_num},{fourth_num}\n")
            elif ',same' in line and ',same' in next_line:
                nums = line.split(',')
                third_num = nums[2]
                fourth_num = nums[3]
                output_lines.append(line)
                output_lines.append(f"{third_num},{fourth_num}\n")
            elif ',same' in line and ',different' in next_line:
                nums = line.split(',')
                first_num = nums[0]
                last_num = nums[5]
                last_num_before_zero = last_num.split('0')[0]
                output_lines.append(line)
            elif ',different' in line and ',different' in next_line:
                nums = line.split(',')
                first_num = nums[0]
                if len(nums) > 5:
                    last_num = nums[5]
                else:
                    last_num = nums[1]
                last_num_before_zero = last_num.split('0')[0]
                output_lines.append(line)
                output_lines.append(f"{first_num},{last_num}\n")
            else:
                output_lines.append(line)
        output_lines.append(lines[-1])
        output = os.path.join(chain_dir,"output.txt")
        with open(output, 'w') as f:
            f.writelines(output_lines)
        ##################################################################################################
        output1 = os.path.join(chain_dir, "output1.txt")
        with open(output, 'r') as f_in, open(output1, 'w') as f_out:
            lines = f_in.readlines()
            i = 0
            line = lines[i]
            nums = line.split(',')
            while i+2 < len(lines):
                line = lines[i]
                second_line = lines[i+1] if i+1 < len(lines) else None
                third_line = lines[i+2] if i+2 < len(lines) else None
                if ',different' in line and ',same' not in second_line and ',different' in third_line :
                    line_list = line.strip().split(',')
                    second_line_list = second_line.strip().split(',') if second_line else None
                    third_line_list = third_line.strip().split(',') if third_line else None
                    if len(line_list) > 5:
                        del line_list[1:5]
                        second_line_list = third_line.strip().split(',') if third_line else second_line_list
                        f_out.write(','.join(line_list)+'\n')
                        f_out.write(','.join(second_line_list)+'\n')
                        i += 3
                    elif len(line_list) < 5:
                        del second_line_list[0:5]
                        f_out.write(','.join(line_list)+'\n')
                        i += 2
                else:
                    f_out.write(line)
                    i += 1
            if i < len(lines):
                f_out.write(lines[i])
            if i+1 < len(lines):
                f_out.write(lines[i+1])
        ###################################################################################################
        output2 = os.path.join(chain_dir, "output2.txt")
        with open(output1, 'r') as f_in, open(output2, 'w') as f_out:
            lines = f_in.readlines()
            i = 0
            line = lines[i]
            nums = line.split(',')
            while i+2 < len(lines):
                line = lines[i]
                second_line = lines[i+1] if i+1 < len(lines) else None
                third_line = lines[i+2] if i+2 < len(lines) else None
                if ',different' in line and ',same' not in second_line and ',different' in third_line :
                    line_list = line.strip().split(',')
                    second_line_list = second_line.strip().split(',') if second_line else None
                    third_line_list = third_line.strip().split(',') if third_line else None
                    if len(line_list) > 5:
                        del line_list[1:5]
                        second_line_list = third_line.strip().split(',') if third_line else second_line_list                    
                        f_out.write(','.join(line_list)+'\n')
                        f_out.write(','.join(second_line_list)+'\n')
                        i += 3
                    elif len(line_list) < 5:
                        f_out.write(','.join(line_list)+'\n')
                        f_out.write(','.join(second_line_list)+'\n')
                        i += 2
                else:
                    f_out.write(line)
                    i += 1
            if i < len(lines):
                f_out.write(lines[i])
            if i+1 < len(lines):
                f_out.write(lines[i+1])
        ###################################################################################################
        output3 = os.path.join(chain_dir, "output3.txt")
        with open(output2, 'r') as f_in, open(output3, 'w') as f_out:
            lines = f_in.readlines()
            i = 0
            line = lines[i]
            nums = line.split(',')
            while i +2 < len(lines):
                if len(nums) > 5:
                    line = lines[i]
                    second_line = lines[i+1] if i+1 < len(lines) else None
                    third_line = lines[i+2] if i+2 < len(lines) else None
                    if ',different' in line and ',same' in third_line :
                        line_list = line.strip().split(',')
                        second_line_list = second_line.strip().split(',') if second_line else None
                        third_line_list = third_line.strip().split(',') if third_line else None
                        del line_list[2:4]
                        second_line_list.extend(third_line.strip().split(',')[0:]) if third_line else None
                        del second_line_list[4:6]
                        f_out.write(','.join(line_list)+'\n')
                        f_out.write(','.join(second_line_list)+'\n')
                        i += 3
                    else:
                        f_out.write(line)
                        i += 1
            if i < len(lines):
                f_out.write(lines[i])
            if i+1 < len(lines):
                f_out.write(lines[i+1])
        ###################################################################################################
        output4 = os.path.join(chain_dir, "output4.txt")
        with open(output3, 'r') as f_in, open(output4, 'w') as f_out:
            lines = f_in.readlines()
            i = 0
            while i < len(lines):
                if i+2 < len(lines):
                    line = lines[i]
                    second_line = lines[i+1]
                    third_line = lines[i+2]
                    if (',same' in second_line and ',different' in third_line) and (',same' not in line and ',different' not in line):
                        line_list = line.strip().split(',')
                        second_line_list = second_line.strip().split(',')
                        third_line_list = third_line.strip().split(',')
                        line_list.extend(second_line.strip().split(',')[0:])
                        del line_list[3:7]
                        f_out.write(','.join(line_list)+'\n')
                        i += 2
                    else:
                        f_out.write(line)
                        i += 1
                else:
                    f_out.write(lines[i])
                    i += 1
        ###################################################################################################
        output5 = os.path.join(chain_dir, "output5.txt")
        with open(output4, 'r') as f_in, open(output5, 'w') as f_out:
            lines = f_in.readlines()
            i = 0
            while i < len(lines):
                if i+2 < len(lines):
                    line = lines[i]
                    second_line = lines[i+1]
                    third_line = lines[i+2]
                    if ',same' in line and ',different' not in second_line and ',same' in third_line:
                        line_list = line.strip().split(',')
                        second_line_list = second_line.strip().split(',')
                        third_line_list = third_line.strip().split(',')
                        second_line_list.extend(third_line.strip().split(',')[0:])
                        del second_line_list[4:6]
                        f_out.write(','.join(line_list)+'\n')
                        f_out.write(','.join(second_line_list)+'\n')
                        i += 3
                    else:
                        f_out.write(line)
                        i += 1
                else:
                    f_out.write(lines[i])
                    i += 1
        ###################################################################################################
        output6 = os.path.join(chain_dir, "output6.txt")
        with open(output5, 'r') as f_in, open(output6, 'w') as f_out:
            lines = f_in.readlines()
            i = 0
            while i < len(lines):
                if i+2 < len(lines):
                    line = lines[i]
                    second_line = lines[i+1]
                    third_line = lines[i+2]
                    if ',same' in line and (',different' not in second_line and ',same' not in second_line) and ',same' in third_line:
                        line_list = line.strip().split(',')
                        second_line_list = second_line.strip().split(',')
                        third_line_list = third_line.strip().split(',')
                        second_line_list.extend(third_line.strip().split(',')[0:])
                        del second_line_list[4:6]
                        f_out.write(','.join(line_list)+'\n')
                        f_out.write(','.join(second_line_list)+'\n')
                        i += 2
                    else:
                        f_out.write(line)
                    i += 1
                else:
                    f_out.write(lines[i])
                    i += 1
        ###################################################################################################
        output7 = os.path.join(chain_dir, "output7.txt")
        with open(output6, 'r') as f_in, open(output7, 'w') as f_out:
            lines = f_in.readlines()
            i = 0
            while i < len(lines):
                if i+2 < len(lines):
                    line = lines[i]
                    second_line = lines[i+1]
                    third_line = lines[i+2]
                    if (',same' in second_line and (',different' in third_line if len(third_line) > 0 else True)) and (',different' in line):
                        line_list = line.strip().split(',')
                        second_line_list = second_line.strip().split(',')
                        third_line_list = third_line.strip().split(',')
                        del second_line_list[3:5]
                        f_out.write(','.join(line_list)+'\n')
                        f_out.write(','.join(second_line_list)+'\n')
                        i += 2
                    else:
                        f_out.write(line)
                        i += 1
                else:
                    f_out.write(lines[i])
                    i += 1
        ###################################################################################################
        output8 = os.path.join(chain_dir, "output8.txt")
        with open(output7, 'r') as f_in, open(output8, 'w') as f_out:
            lines = f_in.readlines()
            i = 0
            while i + 1 < len(lines):
                line = lines[i]
                second_line = lines[i+1] if i+1 < len(lines) else None
                third_line = lines[i+2] if i+2 < len(lines) else None
                if third_line is None and ',same' in second_line :
                    line_list = line.strip().split(',')
                    second_line_list = second_line.strip().split(',') if second_line else None
                    third_line_list = third_line.strip().split(',') if third_line else None
                    del second_line_list[3:5]
                    f_out.write(','.join(line_list)+'\n')
                    f_out.write(','.join(second_line_list)+'\n')
                    i += 2
                else:
                    f_out.write(line)
                    i += 1
            if i < len(lines):
                f_out.write(lines[i])
            if i+1 < len(lines):
                f_out.write(lines[i+1])
        ###################################################################################################
        output9 = os.path.join(chain_dir, "output9.txt")
        with open(output8, 'r') as f_in, open(output9, 'w') as f_out:
            lines = f_in.readlines()
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                line_split = line.split(',')
                if len(line_split) >= 2:
                    first_two = ','.join(line_split[:2])  
                    next_line = lines[i + 1].strip() if i + 1 < len(lines) else None
                    if next_line and next_line.startswith(first_two):
                        if 'different' in line or 'same' in line:
                            f_out.write(line + '\n')  
                        if 'different' in next_line or 'same' in next_line:
                            f_out.write(next_line + '\n')  
                        i += 2  
                    else:
                        f_out.write(line + '\n')  
                        i += 1
                else:
                    f_out.write(line + '\n')  
                    i += 1
        ###################################################################################################
        output_file = os.path.join(chain_dir, "output_file.txt")
        with open(output9, 'r') as f_in, open(output_file, 'w') as f_out:
            lines = f_in.readlines()
            lines = [line for line in lines if line.strip()]
            for line in lines:
                nums = line.strip().split(',')
                nums = [int(num) for num in nums if num != 'different' and num != 'same']
                n = len(nums)
                f_out.write('{:16d}'.format(nums[0]))
                for i in range(1, n, 6):
                    f_out.write(''.join('{:7d}'.format(num) for num in nums[i:i+6])+'\n')
        ###################################################################################################
        output_file1 = os.path.join(chain_dir, "output_file1.txt")
        with open(output_file, 'r') as infile, open(output_file1, 'w') as outfile:
            for line in infile:
                numbers = list(map(int, line.split()))
                if len(numbers) >= 2 and numbers[1] - numbers[0] > 1:
                    numbers[1] = -numbers[1]
                outfile.write('\t'.join(map(str, numbers)) + '\n')
        coordinates = os.path.join(chain_dir, "coordinates.txt")
        with open(output_file1, 'r') as infile, open(coordinates, 'w') as outfile:
            # outfile.write('      INDAT(1)= 0\n')
            for line in infile:
                numbers = list(map(int, line.split()))
                if len(numbers) >= 3:
                    if numbers[-2] - numbers[-3] > 1:
                        numbers[-2] = -numbers[-2]
                    if numbers[-2] - numbers[-1] > 1:
                        numbers[-1] = -numbers[-1]
                elif len(numbers) == 2 and numbers[-1] != 0:
                    if numbers[1] - numbers[0] > 1:
                        numbers[1] = -numbers[1]
                outfile.write(' ' * 10 + ''.join(['{:>7}'.format(str(x)) for x in numbers]) + '\n')
    ###################################################################################################
    def concatenate_cord_files_sorted(self, chain_dir):
        chain_folders = sorted([d for d in os.listdir(chain_dir) if os.path.isdir(os.path.join(chain_dir, d))])
        output_file_name = os.path.join(chain_dir, 'concatenated_cord.txt')
        with open(output_file_name, 'w') as output_file:
            for chain_folder in chain_folders:
                cord_file_path = os.path.join(chain_dir, chain_folder, 'coordinates.txt')
                if os.path.exists(cord_file_path):
                    with open(cord_file_path, 'r') as cord_file:
                        output_file.write(cord_file.read())
    ###################################################################################################
    def generate_output_file(self, header_content, body_content, tail_content):
        with open(self.output_file, 'a') as out_file:
            out_file.write(header_content)
            with open(self.input_pdb_file, 'r') as pdb_file_FMO:
                lines = pdb_file_FMO.readlines()
            ca_indexes = [i for i, line in enumerate(lines) if line[12:16].strip() == "CA"]
            num_lines = 0
            for i, ca_index in enumerate(ca_indexes[:-1]):
                current_residue_number = int(lines[ca_index][22:26].strip())
                next_residue_number = int(lines[ca_indexes[i + 1]][22:26].strip())
                if current_residue_number + 1 == next_residue_number:
                    if (lines[ca_index - 1][12:16].strip() == "N" and
                        lines[ca_index + 1][12:16].strip() == "C" and
                        lines[ca_index + 2][12:16].strip() == "O"):
                        ca_number = int(lines[ca_index][6:11].strip())
                        c_number = int(lines[ca_index + 1][6:11].strip())
                        num_lines += 1
            out_file.write(f"\n      MAXBND={int(num_lines) + 1}\n")
    ###################################################################################################
            with open(self.input_pdb_file, 'r') as pdb_file_FMO:
                res_num = None
                res_name = None
                res_chain_id = None
                res_list = []
                charge_list = []
                res_per_line = 5
                charge_per_line = 10
                line_count = 0
                for line in pdb_file_FMO:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        this_res_chain_id = line[21:22].strip()
                        this_res_num = line[22:26].strip()
                        this_res_name = line[17:20].strip()
                        this_res_charge = line[77:80].strip()
                        if this_res_name not in ['ACE', 'NMA']:
                            if this_res_num.isdigit():
                                res_num_str = str(int(this_res_num))
                            else:
                                res_num_str = this_res_num
                            if this_res_num != res_num or this_res_name != res_name:
                                if res_num is not None and res_name is not None:
                                    res_list.append(res_name + res_chain_id + res_num)
                                    charge_list.append(charge)
                                charge = 0
                            res_num = this_res_num
                            res_name = this_res_name
                            res_chain_id = this_res_chain_id
                            if this_res_charge == 'O1+' or this_res_charge == 'N1+' or this_res_charge == 'n2+' or this_res_charge == 'g2+':
                                charge += 1
                            elif this_res_charge == 'O1-' or this_res_charge == 'S1-' or this_res_charge == 'l1-' or this_res_charge == 'N1-':    
                                charge -= 1
                res_list.append(res_name + res_chain_id + res_num_str)
                charge_list.append(charge)
                for i in range(len(charge_list)):
                    if charge_list[i] == -2:
                        charge_list[i] = -1 
                    elif charge_list[i] == 2:
                        charge_list[i] = 1 
                lig_info = self.input_pdb_file.split('_')[-6:-3]
                if lig_info[-1].isdigit() and len(lig_info[-1]) == 1:
                    search_string = f"{lig_info[0]} {lig_info[1]}   {lig_info[2]}"
                elif lig_info[-1].isdigit() and len(lig_info[-1]) == 2:
                    search_string = f"{lig_info[0]} {lig_info[1]}  {lig_info[2]}"
                elif lig_info[-1].isdigit() and len(lig_info[-1]) == 3:
                    search_string = f"{lig_info[0]} {lig_info[1]} {lig_info[2]}"
                else:
                    search_string = ' '.join(lig_info)
                LIG_charge = None
                with open('../FMOPhore.log', 'r') as log_file:
                    for line in log_file:
                        if search_string in line:
                            parts = line.split()
                            LIG_charge = int(parts[-1])  
                if LIG_charge is not None:
                    charge_list[-1] = LIG_charge  
                nfrag = len(res_list)
                out_file.write('''      NLAYER=1\n''')
                out_file.write(f'      NFRAG={nfrag}\n')
                out_file.write('      ICHARG(1)= ')
                for i, charge in enumerate(charge_list):
                    if i > 0 and i % charge_per_line == 0:  
                        out_file.write('\n                 ')
                    if i == len(charge_list) - 1:  
                        out_file.write(f'{charge:2d}')
                    else:
                        out_file.write(f'{charge:2d}, ')
                out_file.write('\n')
                out_file.write('      FRGNAM(1)= ')
                for i, res in enumerate(res_list):
                    res_name, res_chain_id, res_num = res[:3], res[3], res[4:]
                    if res_name == 'LIG':
                        res = res_name
                    if i > 0 and i % res_per_line == 0:
                        out_file.write('\n                 ')
                    if i == len(res_list) - 1:
                        out_file.write(res)
                    else:
                        out_file.write(res + ', ')
                out_file.write('\n')
                # print(f"Output written to {self.output_file}")
        ###################################################################################################
            with open('chains/concatenated_cord.txt', 'r') as f:
                coordinates = f.read().rstrip()
                out_file.write('      INDAT(1)= 0\n')
                out_file.write(coordinates)
            shutil.rmtree("chains")
        ###################################################################################################
            with open(self.input_pdb_file, 'r') as pdb_file_FMO:
                lines = pdb_file_FMO.readlines()
            ligand_first, ligand_last = None, None
            for line in lines:
                if line.startswith("HETATM"):    
                    atom_number = int(line[6:11].strip())
                    if ligand_first is None or atom_number < ligand_first:
                        ligand_first = int(atom_number)
                    if ligand_last is None or atom_number > ligand_last:
                        ligand_last = int(atom_number)
            if ligand_first is not None and ligand_last is not None:
                out_file.write(f"\n{ligand_first:>17} {-ligand_last:>6}      0\n")
            else:
                print("Could not find ligand in the PDB file")
        ###################################################################################################
            out_file.write(body_content + "\n")
        ###################################################################################################
            with open(self.input_pdb_file, 'r') as pdb_file_FMO:
                lines = pdb_file_FMO.readlines()
            ca_indexes = [i for i, line in enumerate(lines) if line[12:16].strip() == "CA"]
            num_lines = 0
            for i, ca_index in enumerate(ca_indexes[:-1]):
                current_residue_number = int(lines[ca_index][22:26].strip())
                next_residue_number = int(lines[ca_indexes[i+1]][22:26].strip())    
                if current_residue_number + 1 == next_residue_number:
                    if (lines[ca_index-1][12:16].strip() == "N" and
                        lines[ca_index+1][12:16].strip() == "C" and
                        lines[ca_index+2][12:16].strip() == "O"):
                        ca_number = int(lines[ca_index][6:11].strip())
                        c_number = int(lines[ca_index+1][6:11].strip())
                        out_file.write(f"{'-' + str(ca_number):>10} {c_number:>9}  DFTB-C\n")
                        num_lines += 1
            out_file.write(''' $END
 $DATA\n''')
        ###################################################################################################
            out_file.write(f''' FMO calculation : {self.input_pdb_file}\n''')
            out_file.write(tail_content + "\n")
        ###################################################################################################
            with open(self.input_pdb_file, 'r') as f:
                lines = f.readlines()
            out_lines = []
            for line in lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_num = line[6:11].strip()
                    if line.startswith('HETATM'):
                        atom_num = str(int(atom_num)) 
                    atom_name = line[12:16].strip()
                    if atom_name[:2].upper() in ['CL', 'BR', 'MG']:
                        atom_name = atom_name[:2].upper()
                    elif atom_name[0].isdigit():
                        atom_name = 'H' 
                    else:
                        atom_name = atom_name[0].upper()  
                    x = line[30:38].strip()
                    y = line[38:46].strip()
                    z = line[46:54].strip()
                    out_lines.append(f'{atom_num:>7} {atom_name:>6} {x:>16} {y:>12} {z:>12}\n')
            out_file.write(''.join(out_lines))
            out_file.write(' $END\n') 
###################################################################################################
if __name__ == "__main__":
    """
    FMOPhore V.0.1 - Protein-ligand complex Fragmention - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    pdb_processor = Fragmention_Processor(self.pdb_file, self.qm)
    pdb_processor.process_Fragmention()