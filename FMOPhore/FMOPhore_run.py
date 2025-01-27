import sys
import os
import argparse
import subprocess
import shutil, os, glob
import time
import re
import threading
from multiprocessing import Pool
import multiprocessing
import getpass
import timeout_decorator
from datetime import datetime
from tqdm import tqdm
import concurrent.futures

from FMOPhore import merge_ligs_prot, download_pdb, download_pdbs_from_file
from FMOPhore import PDBProcessor
from FMOPhore import aligner
from FMOPhore import LIGProcessor
from FMOPhore import CutoffProcessor
from FMOPhore import SplitProcessor
from FMOPhore import Pharmacophores
from FMOPhore import Analysis
from FMOPhore import Fragmention_Processor
from FMOPhore import ComplexAnalysis
from FMOPhore import LibAnalysis
####################################################################################################
def parse_cutoff(value):
    if value == 'no_cutoff':
        return sys.maxsize
    else:
        return int(value)
####################################################################################################
# Parse command-line arguments
####################################################################################################
mandatory_params_text = """
\033[1mMandatory parameters:\033[0m

    -dir,       --directory             : Process all PDB files in a directory
    -com,       --Prot_complex          : Process a single complex PDB file .pdb
    -PDB,       --PDB_ID                : PDB ID or .txt file containing PDB IDs.
                                          Process a single complex PDB directly from Protein Data Bank > https://www.rcsb.org/ 
    -prot       --protein_pdb_file      : Path to protein PDB file
    -ligs       --ligand_files          : Path to single ligand PDB file or directory of ligand PDB files

    -qm         --qm_calculation        : MP2 or DFTB (GAMESS software required/ https://www.msg.chem.iastate.edu/gamess/download.html)
    -d          --distance_cutoff       : Distance cutoff for selecting residues. Use "no_cutoff" to select the whole protein.
    
    -t          --timer                 : Timer in days e.g.: 1 day ==> -t 1  
    -c          --cpus                  : Please specify the number of cpus (cores; number of jobs to parallelize)
"""
optional_params_text = """
\033[1mOptinal parameters:\033[0m

    -DA-FMO,    --trajectory            : To do DA-FMO, read trajectory "traj.dcd" 
                                        : Dy-FMOPhore requires SuMD (https://github.com/molecularmodelingsection/SuMD) 
                                        :                      ACEMD (https://software.acellera.com/acemd/tutorial.html)
    -cof,       --cofactor              : If you have co-factor in the system, provide the name in this form, e.g.: LYS-600'
    -BE,        --Binding_Energy        : Calculate the Binding Energy (ΔE) = (E_Complex) − (E_Protein) − (E_Ligand)
    -lib,       --same_target           : If same target and different ligands
    -align,     --align                 : If needed to align structures
    -analysis   --FMOPhore_analysis     : To run analysis only if calculations has been completed.
"""
parser = argparse.ArgumentParser(description=':||: FMOPhore 0.1 :||:')
mandatory_group = parser.add_argument_group(mandatory_params_text)
mandatory_group.add_argument('-These are required parameters', '-!', type=str, help='')
mandatory_group.add_argument('-dir', '--directory', type=str, help=argparse.SUPPRESS)
mandatory_group.add_argument('-com', '--Prot_complex', type=str, help=argparse.SUPPRESS)
mandatory_group.add_argument('-PDB', '--PDB_ID', type=str, default=None, help=argparse.SUPPRESS)
mandatory_group.add_argument('-prot', '--protein_pdb_file', type=str, help=argparse.SUPPRESS)
mandatory_group.add_argument('-ligs', '--ligand_files', type=str, help=argparse.SUPPRESS)
mandatory_group.add_argument('-qm', '--qm_calculation', type=str, default=None, help=argparse.SUPPRESS)
mandatory_group.add_argument('-d', '--distance_cutoff', type=parse_cutoff, default=None, help=argparse.SUPPRESS, const=sys.maxsize, nargs='?')
mandatory_group.add_argument('-t', '--timer', type=int, default=None, help=argparse.SUPPRESS)
mandatory_group.add_argument('-c', '--cpus', type=int, default=None, help=argparse.SUPPRESS)
optional_group = parser.add_argument_group(optional_params_text)
optional_group.add_argument('-These are optinal parameters', '-...', type=str, help='')
optional_group.add_argument('-DA-FMO', '--trajectory', type=str, default=None, help=argparse.SUPPRESS)
optional_group.add_argument('-cof', '--cofactor', type=str, default=None, help=argparse.SUPPRESS)
optional_group.add_argument('-BE', '--Binding_Energy', action='store_true', default=None, help=argparse.SUPPRESS)
optional_group.add_argument('-lib', '--same_target', action='store_true', default=None, help=argparse.SUPPRESS)
optional_group.add_argument('-align', '--align', action='store_true', default=None, help=argparse.SUPPRESS)
optional_group.add_argument('-analysis', '--FMOPhore_analysis', action='store_true', default=None, help=argparse.SUPPRESS)
parser.epilog = "Peter E.G.F. Ibrahim."
args = parser.parse_args() 

if args.qm_calculation:
    print('''Please install GAMESS software: 
        https://www.msg.chem.iastate.edu/gamess/download.html
        If you have it installed and running:
        Recommended a GPU cluster equipped with at least 2 GPUs and 20 CPU.
        Contact the author Peter E.G.F. Ibrahim: 2448959@dundee.ac.uk, for further details on running FMOPhore to its full potential.''')


if (args.directory or args.Prot_complex or args.PDB_ID or (args.protein_pdb_file and args.ligand_files)):
    pass
else:
    parser.print_help()
    sys.exit(1)
if args.distance_cutoff == 'no_cutoff':
    distance_cutoff = None
else:
    distance_cutoff = args.distance_cutoff
if distance_cutoff is None:
    distance_cutoff = sys.maxsize
if args.cofactor:
    try:
        cofactor_resname, cofactor_resnum = args.cofactor.split('-')
    except ValueError:
        raise ValueError('Invalid cofactor argument: %s' % args.cofactor)
else:
    cofactor_resname = None
    cofactor_resnum = None
script_dir = os.path.dirname(os.path.abspath(__file__))
####################################################################################################
# change the name of files if they have a . to _
####################################################################################################
folder_path = './'
pdb_files = glob.glob(os.path.join(folder_path, '*.pdb'))
for pdb_file in pdb_files:
    file_name     = os.path.splitext(os.path.basename(pdb_file))[0]
    new_file_name = file_name.replace('.', '_') + '.pdb'
    new_file_path = os.path.join(folder_path, new_file_name)
    os.rename(pdb_file, new_file_path)
####################################################################################################
# Arguments
####################################################################################################
def copy_and_update_rungms(target_directory, files_to_copy):
    package_path = os.path.dirname(__file__)  
    for file_name in files_to_copy:
        source = os.path.join(package_path, "QM_run", file_name)
        destination = os.path.join(target_directory, file_name)
        if not os.path.exists(source):
            raise FileNotFoundError(f"Source file does not exist: {source}")
        shutil.copy(source, destination)
        if file_name in ["DFTB", "MP2"]:  
            with open(destination, "r") as f:
                content = f.read()
            updated_content = content.replace(
                "/FMOPhore/QM_runs/", f"{package_path}/QM_run/")
            with open(destination, "w") as f:
                f.write(updated_content)

if args.qm_calculation == "DFTB":
    files_to_copy = ["run_DFTB.py", "rungms", "DFTB"]
elif args.qm_calculation == "MP2": 
    files_to_copy = ["run_MP2.py", "rungms", "MP2"]
####################################################################################################
############################ To define the PDBs that FMOPhore will work on. #########################
pdb_files = []

if args.Prot_complex:
    ref_file  = args.Prot_complex
    pdb_files = [ref_file]
if args.directory and not (args.protein_pdb_file and args.ligand_files):
    pdb_files = glob.glob(os.path.join("*.pdb"))
if args.PDB_ID:
    error_pdb_ids = []
    if args.PDB_ID.endswith('.txt'):
        print(f"Processing PDB list from file: {args.PDB_ID}")
        download_pdbs_from_file(args.PDB_ID, error_pdb_ids)
    else:
        print(f"Downloading single PDB file: {args.PDB_ID}")
        download_pdb(args.PDB_ID, error_pdb_ids)
    pdb_files = glob.glob("*.pdb")
    if error_pdb_ids:
        with open("error_pdb_ids.txt", "w") as error_file:
            error_file.write("\n".join(error_pdb_ids))
        print(f"Failed downloads recorded in error_pdb_ids.txt")
    print(f"Downloaded PDB files: {pdb_files}")
if args.protein_pdb_file and args.ligand_files and not args.directory:
    print("Merging protein and ligand PDB files...")
    output_dir = "./"
    if os.path.isfile(args.ligand_files):
        merge_ligs_prot(args.protein_pdb_file, args.ligand_files, output_dir)
    elif os.path.isdir(args.ligand_files):
        ligand_files = glob.glob(os.path.join(args.ligand_files, "*.pdb"))
        for ligand_file in ligand_files:
            merge_ligs_prot(args.protein_pdb_file, ligand_file, output_dir)
    pdb_files = glob.glob(os.path.join(output_dir, "*.pdb"))
###################################################################################################
# Library
###################################################################################################
if args.same_target:
    new_directory = "library_analysis/residues"
    os.makedirs(new_directory, exist_ok=True)
def lib_pdbs(pdb_file):
    file_list = glob.glob('library_analysis/residues/*_residues.txt')
    content_list = []
    # Iterate over each file
    for file_name in file_list:
        with open(file_name, 'r') as file:
            content = file.readlines()
            content_list.extend(content)
    content_list = [line.strip().split('-') for line in content_list]
    content_list = [(line[0], line[1], int(line[2])) for line in content_list]
    def custom_sort_key(residue):
        return (residue[1], residue[2])
    sorted_content = sorted(content_list, key=custom_sort_key)
    unique_content = []
    seen = set()
    for line in sorted_content:
        line_str = f"{line[0]}-{line[1]}-{line[2]}"
        if line_str not in seen:
            unique_content.append(line_str)
            seen.add(line_str)
    with open('library_analysis/residues/residues.txt', 'w') as merged_file:
        for line in unique_content:
            merged_file.write(line + '\n')
######################################### Command runner ###########################################
def run_command(command):
    return subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
######################################### PDBProcessor #############################################
excluded_strings = [ 'EOH', 'GZ6','DMS', 'FMT', 
                     'ANISOU', 'PEG', 'GOL', 'BO3', 'EDO', 'BU3',
                     'SO4', 'ACE', 'NMA', 'NME', 'ACT', 'MES', 'MLA',
                     'OCY', 'SEP', 'TPO', 'IPA', 'TRS', 'CO', 'NO3', 'NAG' , 'PO4',
                     'BME', 'CSO', 'IMD', 'TCE', 'MDP', 'IOD', 'PTR','NI', 'TRS'
                     , 'SAH', 'MTA', 'DMS'
                     ]
####################################################################################################
def preparePDBs(pdb_file):
    folder_name = os.path.splitext(pdb_file)[0]
    os.makedirs(folder_name, exist_ok=True)
    # shutil.copy(pdb_file, folder_name)
    shutil.move(pdb_file, folder_name)
    FMOPhore_log = "FMOPhore.log"
    FMOPhore_input_file = os.path.join(script_dir, FMOPhore_log)
    shutil.copy(FMOPhore_input_file, folder_name)
    os.chdir(folder_name)    
    ############################################################
    pdb_processor = PDBProcessor(pdb_file)
    pdb_processor.process_complex()
    os.chdir("..")
    #############################################################
def preparePDBs_all(files):
    with Pool(processes=args.cpus) as pool:    
        pool.map(preparePDBs, files)
######################################### Align library ############################################
def align_PDBs(pdb_info):
    pdb_file, ref_file = pdb_info
    ref_folder_name = os.path.dirname(ref_file)
    parent_dir = os.getcwd()
    ref_base_name = os.path.splitext(os.path.basename(ref_file))[0]
    ref_file_modified = f"{ref_base_name}.pdb"
    def align(pdb_file):
        pdb_folder_name = os.path.dirname(pdb_file)
        pdb_base_name = os.path.splitext(os.path.basename(pdb_file))[0]
        target_dir = os.path.join(pdb_folder_name, pdb_base_name)
        if os.path.exists(target_dir):
            os.chdir(target_dir)
            ref_path = os.path.join(parent_dir, ref_base_name, f"{ref_base_name}.pdb")
            aligner(
                reference_pdb_file=ref_path,
                pdb_files_pattern=f"{pdb_base_name}.pdb",
                output_dir="align"
            )
        else:
            print(f"Directory {target_dir} does not exist to align")
    align(pdb_file)
    align_dir = "align"
    if os.path.exists(align_dir):
        shutil.rmtree(align_dir)
    os.chdir("..")
    #############################################################
def align_library(ref_file, pdb_files_notref):
    tasks = [(pdb_file, ref_file) for pdb_file in pdb_files_notref]
    with Pool(processes=args.cpus) as pool:
        pool.map(align_PDBs, tasks)
######################################### Lig_separate #############################################
def sep_LIGsPDBs(pdb_file):
    folder_name = os.path.splitext(pdb_file)[0]
    os.chdir(folder_name)
    pdb_processor = LIGProcessor(pdb_file)
    ligands = pdb_processor.process_complex_LIG()
    os.chdir("..")
    return (pdb_file, ligands)

def sep_LIGsPDBs_all(files):
    with Pool(processes=args.cpus) as pool:
        file_ligands = pool.map(sep_LIGsPDBs, files)
######################################### Cutoff - distance ########################################
def cutoffPDBs(pdb_file):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    folder_name = os.path.splitext(pdb_file)[0]
    # full_path = os.path.join(os.getcwd(), base_name)
    # print(f"Base name: {base_name}")
    # print(f"Full path: {full_path}")
    os.chdir(folder_name)
    ############################################################
    def cutoffPDB_processor(filename, distance, same_target=False):
        if not any(exclude_str in f"{ligand_info}" for exclude_str in excluded_strings):
            with os.scandir(f"{filename}"):
                os.chdir(filename)
                processor = CutoffProcessor(
                    pdb_file=f"{filename}.pdb",
                    distance_cutoff=args.distance_cutoff,
                    same_target=args.same_target,
                )
                processor.cutoff_complex()
                print(f"cutoffPDB_processor {filename}.pdb -d {distance}")
                os.chdir('..')
    ligands_by_chain = {}
    with open("FMOPhore.log", 'r') as f_in:
        for line in f_in:
            line = line.strip()
            if "falls in chain" in line:
                ligand_info = line[0:9]
                ligand_info = "_".join(ligand_info.split())
                chain_id = line[4:5]  
                if chain_id not in ligands_by_chain:
                    ligands_by_chain[chain_id] = []
                ligands_by_chain[chain_id].append(f"{ligand_info}")
    if any(ligands for ligands in ligands_by_chain.values()):
        for chain_id, ligands in ligands_by_chain.items():
            for ligand in ligands:
                ligand_info = "_".join(ligand.split())
                target_dir = f"{base_name}_{ligand_info}"
                if os.path.exists(target_dir):
                    cutoffPDB_processor(target_dir, args.distance_cutoff, args.same_target)
                else:
                    print(f"Directory {target_dir} does not exist to cutoff")
    os.chdir("../")
    #############################################################
def cutoff_PDBs_all(files):
    with Pool(processes=args.cpus) as pool:    
        pool.map(cutoffPDBs, files)
######################################### SPLIT QM_FMOPhore ########################################
def splitPDBs(pdb_file):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    folder_name = os.path.splitext(pdb_file)[0]
    os.chdir(folder_name)
    ############################################################
    def splitProcessor(filename, distance, same_target=False):
        if not any(exclude_str in f"{ligand_info}" for exclude_str in excluded_strings):
            with os.scandir(f"{filename}"):
                os.chdir(f"{filename}")
                processor = SplitProcessor(
                    pdb_file=f"{filename}.pdb",
                    distance_cutoff=args.distance_cutoff,
                    same_target=args.same_target,
                    Binding_Energy=args.Binding_Energy)
                processor.split_complex()
                print(f"splitProcessor {filename}.pdb")
                os.chdir('..')
    ligands_by_chain = {}
    with open("FMOPhore.log", 'r') as f_in:
        for line in f_in:
            line = line.strip()
            if "falls in chain" in line:
                ligand_info = line[0:9]
                ligand_info = "_".join(ligand_info.split())
                chain_id = line[4:5]  
                if chain_id not in ligands_by_chain:
                    ligands_by_chain[chain_id] = []
                ligands_by_chain[chain_id].append(f"{ligand_info}")
    if any(ligands for ligands in ligands_by_chain.values()):
        for chain_id, ligands in ligands_by_chain.items():
            for ligand in ligands:
                ligand_info = "_".join(ligand.split())
                target_dir = f"{base_name}_{ligand_info}"
                if os.path.exists(target_dir):
                    splitProcessor(target_dir, args.distance_cutoff, args.same_target)
                else:
                    print(f"Directory {target_dir} does not exist")
    os.chdir("../")
    #############################################################
def splitPDBs_all(files):
    with Pool(processes=args.cpus) as pool:    
        pool.map(splitPDBs, files)
######################################### SPLIT QM_FMOPhore ########################################
def FMOPhore(pdb_file):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    folder_name = os.path.splitext(pdb_file)[0]
    os.chdir(folder_name)
    #############################################################
    def run_FMOPhore(filename, ligand_info, distance):
        if not any(exclude_str in f"{ligand_info}" for exclude_str in excluded_strings):
            with os.scandir(f"{os.path.basename(pdb_file[:-4])}_{ligand_info}"):
                os.chdir(f"{os.path.basename(pdb_file[:-4])}_{ligand_info}")
                processor = Pharmacophores(
                    pdb_file=f"{os.path.basename(pdb_file[:-4])}_{ligand_info}.pdb")
                processor.run_Pharmacophores()
                os.chdir('..')
                print(os.getcwd())
        os.chdir('..')
    ligands_by_chain = {}
    with open("FMOPhore.log", 'r') as f_in:
        for line in f_in:
            line = line.strip()
            if "falls in chain" in line:
                ligand_info = line[0:9]
                ligand_info = "_".join(ligand_info.split())
                chain_id = line[4:5]  
                if chain_id not in ligands_by_chain:
                    ligands_by_chain[chain_id] = []
                ligands_by_chain[chain_id].append(f"{ligand_info}")
    if any(ligands for ligands in ligands_by_chain.values()):
        for chain_id, ligands in ligands_by_chain.items():
            for ligand in ligands:
                ligand_info = "_".join(ligand.split())
                if not any(exclude_str in f"{ligand_info}" for exclude_str in excluded_strings):
                    target_dir = f"{base_name}_{ligand_info}"
                    if os.path.exists(target_dir):
                        run_FMOPhore(f"{os.path.basename(pdb_file[:-4])}_{ligand_info}", ligand_info, f"{args.distance_cutoff}")
                    else:
                        print(f"Directory {target_dir} does not exist")
    os.chdir("../")
    #############################################################
def FMOPhore_all(files):
    with Pool(processes=args.cpus) as pool:    
        pool.map(FMOPhore, files)
######################################### QM_FMOPhore ##############################################
def QM_run(pdb_file):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    folder_name = os.path.splitext(pdb_file)[0]
    os.chdir(folder_name)  
    def run_QM(filename, ligand_info, distance):
        if not any(exclude_str in ligand_info for exclude_str in excluded_strings):
            ligand_dir = f"{os.path.basename(pdb_file[:-4])}_{ligand_info}"
            if os.path.exists(ligand_dir):
                os.chdir(ligand_dir) 
                for file_to_copy in files_to_copy:
                    source_file = os.path.join(script_dir, "QM_run", file_to_copy)
                    destination_file = os.path.join(file_to_copy)
                    shutil.copy(source_file, destination_file)
                    target_directory = os.getcwd()
                def wait_for_file(file_path, interval=20):
                    while not os.path.exists(file_path):
                        time.sleep(interval)
                def process_qm_calculation(pdb_file, ligand_info, distance, calculation_type, processor_class):
                    output_file = f"{os.path.basename(pdb_file[:-4])}.inp"
                    processor = processor_class(input_pdb_file=pdb_file, output_file=output_file, qm=calculation_type)
                    processor.process_Fragmention()
                def process_binding_energy(pdb_file, ligand_info, calculation_type, prot_processor_class, lig_processor_class):
                    protein_file = f"{os.path.basename(pdb_file[:-4])}.inp"
                    processor = prot_processor_class(input_pdb_file=pdb_file, output_file=protein_file, qm=calculation_type)
                    processor.process_prot_Fragmention()
                    lig_file = f"{os.path.basename(pdb_file[:-4])}.inp"
                    processor = lig_processor_class(input_pdb_file=pdb_file, output_file=lig_file, qm=calculation_type)
                    processor.process_lig_Fragmention()
                com_pdb_file = f"{os.path.basename(pdb_file[:-4])}_{ligand_info}_FMO_{distance}_H.pdb"
                wait_for_file(com_pdb_file)
                if args.qm_calculation == "DFTB":
                    copy_and_update_rungms(target_directory, files_to_copy)
                    process_qm_calculation(com_pdb_file, ligand_info, distance, "DFTB", Fragmention_Processor)
                elif args.qm_calculation == "MP2":
                    copy_and_update_rungms(target_directory, files_to_copy)
                    process_qm_calculation(com_pdb_file, ligand_info, distance, "MP2", Fragmention_Processor)
                # if args.Binding_Energy:
                #     prot_pdb_file = f"{os.path.basename(pdb_file[:-4])}_{ligand_info}_FMO_protein.pdb"
                #     lig_pdb_file = f"{os.path.basename(pdb_file[:-4])}_{ligand_info}_FMO_lig_H.pdb"
                #     if args.qm_calculation == "DFTB":
                #         process_binding_energy(com_pdb_file, ligand_info, "DFTB", prot_Fragmention_Processor, lig_Fragmention_Processor)
                #     elif args.qm_calculation == "MP2":
                #         process_binding_energy(com_pdb_file, ligand_info, "MP2", prot_Fragmention_Processor, lig_Fragmention_Processor)
                input_files = [f"{os.path.basename(pdb_file[:-4])}_{ligand_info}_FMO_{args.distance_cutoff}_H.inp"]
                while any(not os.path.exists(file) for file in input_files):
                    time.sleep(20) 
                if args.qm_calculation == "DFTB":
                    command5 = os.path.join("python run_DFTB.py")
                    run_command(command5)
                elif args.qm_calculation == "MP2": 
                    command5 = os.path.join("python run_MP2.py")
                    run_command(command5)

                os.remove('output.pdb')
                os.remove('atom_modified.pdb')
                os.remove(f"{os.path.basename(pdb_file[:-4])}_{ligand_info}_FMO_{args.distance_cutoff}.pdb")    
                os.chdir("..") 
            else:
                print(f"Ligand directory {ligand_dir} does not exist")
    #############################################################
    ligands_by_chain = {}
    with open("FMOPhore.log", 'r') as f_in:
        for line in f_in:
            line = line.strip()
            if "falls in chain" in line:
                ligand_info = "_".join(line[0:9].split())
                chain_id = line[4:5]
                if chain_id not in ligands_by_chain:
                    ligands_by_chain[chain_id] = []
                ligands_by_chain[chain_id].append(ligand_info)
    if any(ligands for ligands in ligands_by_chain.values()):
        for chain_id, ligands in ligands_by_chain.items():
            for ligand in ligands:
                ligand_info = "_".join(ligand.split())
                if not any(exclude_str in ligand_info for exclude_str in excluded_strings):
                    target_dir = f"{base_name}_{ligand_info}"
                    if os.path.exists(target_dir):
                        run_QM(target_dir, ligand_info, args.distance_cutoff)
                    else:
                        print(f"Directory {target_dir} does not exist")
    os.chdir("../")  
###################################################################################################
def MP2(file):
    QM_run(file)
def DFTB(file):
    QM_run(file)
def QM(pdb_files, qm_calculation):
    if args.qm_calculation == "MP2":
        QM_cal = MP2
    elif args.qm_calculation == "DFTB":
        QM_cal = DFTB
    else:
        raise ValueError("Invalid QM calculation type: " + str(args.qm_calculation))
    with Pool(processes=args.cpus) as pool:
        pool.map(QM_cal, pdb_files)
####################################################################################################
################################### QM_FMOPhore analysis ###########################################
####################################################################################################
def QM_FMOPhore_analysis(pdb_file):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    folder_name = os.path.splitext(pdb_file)[0]
    os.chdir(folder_name)
    #############################################################
    def clean_up(files):
        for file_name in files:
            try:
                os.remove(file_name)
            except FileNotFoundError:
                print(f"File {file_name} not found, skipping.")
            except Exception as e:
                pass
    #############################################################
    def analyze_logs(log_files):
        exits = [
        # 'Aborting due',
        'exited gracefully'
         ]
        for log_file in log_files:
            with open(log_file, 'r') as f:
                content = f.read()
                if all(exit_phrase not in content for exit_phrase in exits):
                    return False
        return True
    #############################################################
    #  Analysis
    #############################################################
    def analysis_FMOPhore(filename, 
                          ligand_info, 
                          distance, 
                          same_target=False):
        if not any(exclude_str in f"{ligand_info}" for exclude_str in excluded_strings):
            with os.scandir(f"{filename}"):
                os.chdir(f"{filename}")
                loop_counter = 0
                while loop_counter < 15:
                    log_files = []
                    for dirpath, dirnames, filenames in os.walk('.'):
                        for f in filenames:
                            if f.endswith('.log') and not (f.endswith("FMOPhore.log") or f.endswith("FMOPhore_org.log")):
                                log_files.append(os.path.join(dirpath, f))
                    if len(log_files) == 0:
                        continue
                    all_logs_found = True
                    error_list = [
                        "impossible",
                        " **** ERROR, NO  $DATA   GROUP WAS FOUND", 
                        " DDI Process 0: semget return an error.", 
                        " DDI Process 0: error code 911", 
                        " DDI Process 1: error code 911",
                        " DDI Process 2: error code 911",
                        " DDI Process 3: error code 911",
                        " ddikick.x: Fatal error detected."
                    ]
                    charge_error = " multiplicity"
                    bad_indat = "Bad indat, atom  "
                    for error in error_list:
                        for log_file in log_files:
                            with open(log_file, 'r') as f:
                                file_contents = f.read()
                                if error.strip() in file_contents:
                                    loop_counter += 1
                                    os.chdir(f"{os.path.dirname(log_file)}")
                                    if charge_error.strip() in file_contents:
                                        matches = re.finditer(r"Fragment\s+(\d+):", file_contents)
                                        for match in matches:
                                            fragment_number = int(match.group(1)) - 1 
                                            for inp_file in os.listdir('.'):
                                                if inp_file.endswith('.inp'):
                                                    with open(inp_file, 'r+') as file:
                                                        contents = file.read()
                                                        icharg_pattern = re.compile(r"ICHARG\(1\)=((?:\s*-?\d+,?)+)")
                                                        icharg_match = icharg_pattern.search(contents)
                                                        if icharg_match:
                                                            numbers = re.findall(r"-?\d+", icharg_match.group(1))
                                                            if 0 <= fragment_number < len(numbers):
                                                                if numbers[fragment_number] == '0':
                                                                    numbers[fragment_number] = '-1'
                                                                elif numbers[fragment_number] == '-1':
                                                                    numbers[fragment_number] = '0'
                                                                elif numbers[fragment_number] == '1':
                                                                    numbers[fragment_number] = '0'
                                                                updated_numbers = ', '.join(numbers)
                                                                updated_line = f"ICHARG(1)= {updated_numbers}\n"
                                                                updated_contents = contents[:icharg_match.start()] + updated_line + contents[icharg_match.end():]
                                                                file.seek(0)
                                                                file.write(updated_contents)
                                                                file.truncate()
                                                                # print(f"Updated fragment number {fragment_number + 1} in {inp_file}")                            
                                        def reformat_icharg_line(input_line):
                                            numbers = input_line.split('=')[1].strip().split(',')
                                            numbers = [num.strip() for num in numbers]
                                            formatted_lines = []
                                            for i in range(0, len(numbers), 10):
                                                line_segment = numbers[i:i+10]
                                                formatted_line = ', '.join(f"{int(num):2d}" for num in line_segment)
                                                formatted_lines.append(formatted_line)
                                            output = "      ICHARG(1)=    " + ',\n                    '.join(formatted_lines)
                                            return output
                                        def update_file(filename):
                                            new_content = []
                                            icharg_found = False
                                            with open(filename, 'r') as file:
                                                for line in file:
                                                    if line.strip().startswith("ICHARG(1)="):
                                                        new_content.append(reformat_icharg_line(line))
                                                        icharg_found = True
                                                    else:
                                                        new_content.append(line.rstrip())
                                            if icharg_found:
                                                with open(filename, 'w') as file:
                                                    for line in new_content:
                                                        file.write(line + '\n')
                                        
                                        for dirpath, dirname, filenames in os.walk('.'):
                                            for f in filenames:
                                                if f.endswith('.inp'):
                                                    update_file(f) 
                                    elif bad_indat.strip() in file_contents:
                                        message_FraGrow = bad_indat
                                        with open("../../FMOPhore.log", "a") as log_file:
                                            log_file.write(message_FraGrow)
                                        break              
                                    else:
                                        pass
                                    os.chdir("..")
                                    command6 = f"mv {os.path.dirname(log_file)}/*.inp ./"
                                    command7 = f"rm {os.path.dirname(log_file)}/{os.path.dirname(log_file)}.log"
                                    # result7 = subprocess.run(command7, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                                    run_command(command6)
                                    if args.qm_calculation == "DFTB":
                                        command5 = os.path.join("python run_DFTB.py")
                                        run_command(command5)
                                    elif args.qm_calculation == "MP2": 
                                        command5 = os.path.join("python run_MP2.py")
                                        run_command(command5)
                                    time.sleep(20)   
                                else:
                                    pass      
                #################################################################    
                    log_files = []
                    for dirpath, dirnames, filenames in os.walk('.'):
                        for f in filenames:
                            if f.endswith('.log') and not (f.endswith("FMOPhore.log") or f.endswith("FMOPhore_org.log")):
                                log_files.append(os.path.join(dirpath, f))
                    if len(log_files) == 0:
                        continue
                    all_logs_found = True
                    exits = [
                    # 'Aborting due',
                    'exited gracefully'
                     ]
                    for log_file in log_files:
                        with open(log_file, 'r') as f:
                            content = f.read()
                            if all(exit_phrase not in content for exit_phrase in exits):
                                    all_logs_found = False
                                    break
                    if not log_files:
                        return
                    if analyze_logs(log_files):
                        analysis_dir = "analysis"
                        if not os.path.exists(analysis_dir):
                            os.makedirs(analysis_dir)
                        os.chdir(analysis_dir)
                        run_command(f"cp ../*/*.log ./ && cp ../*/*.inp ./ && cp ../*.pdb ./")
                        
                        pdb_file = f"{filename}_FMO_{distance}_H.pdb"
                        pdb_processor = ComplexAnalysis(pdb_file=pdb_file, same_target=args.same_target)
                        pdb_processor.com_analyze()
                        pdb_processor = Analysis(pdb_file=pdb_file, Binding_Energy=args.Binding_Energy)
                        pdb_processor.analyze()
                        try:
                            print(f"analysis -com ../{pdb_file[:-4]}_{ligand_info}_FMO_{distance}_H.pdb")
                            log_files = glob.glob('*.log')
                            inp_files = glob.glob('*.inp')
                            pdb_files = [f for f in glob.glob('*.pdb') if not (f.startswith('Ring') or f.startswith('linker') or f.startswith('scaffold'))]
                            files_to_clean = log_files + inp_files + pdb_files
                            clean_up(files_to_clean)
                            os.chdir("..")
                            clean_up(['DFTB', 'rungms', 'rungms'])
                            print("parent_dir = " , os.getcwd())
                            os.chdir('..')
                        except OSError as e:
                            if e.errno == 116:
                                print(f"Error changing directory: {e}")
                            else:
                                pass
                        except Exception as e:
                            pass
                        break
    #############################################################
    ligands_by_chain = {}
    with open("FMOPhore.log", 'r') as f_in:
        for line in f_in:
            line = line.strip()
            if "falls in chain" in line:
                ligand_info = line[0:9]
                ligand_info = "_".join(ligand_info.split())
                chain_id = line[4:5]  
                if chain_id not in ligands_by_chain:
                    ligands_by_chain[chain_id] = []
                ligands_by_chain[chain_id].append(f"{ligand_info}")
    if any(ligands for ligands in ligands_by_chain.values()):
        for chain_id, ligands in ligands_by_chain.items():
            for ligand in ligands:
                ligand_info = "_".join(ligand.split())
                if not any(exclude_str in f"{ligand_info}" for exclude_str in excluded_strings):
                    target_dir = f"{base_name}_{ligand_info}"
                    if os.path.exists(target_dir):
                        try:
                            print("analysing:", target_dir)
                            analysis_FMOPhore(f"{os.path.basename(pdb_file[:-4])}_{ligand_info}",
                                                ligand_info, 
                                                f"{args.distance_cutoff}", 
                                                args.same_target)
                        except Exception as e:
                            error_message = f"FMOPhore error occurred: {os.path.basename(pdb_file[:-4])}_{ligand_info}\n"
                            print(error_message)
                            with open("../../../errors.log", "a") as error_file:
                                error_file.write(error_message)  
                            os.chdir("../../")
                        os.chdir("../")
                        print(os.getcwd())
                        return
                        print(f"analysing: Done", target_dir)
                    else:
                        print(f"Directory {target_dir} does not exist")
    os.chdir("..")
    print("parent_dir = " , os.getcwd())
    #############################################################
    message1 = f"FMOPhore - is done on: {folder_name}\n"
    with open("FMOPhore.log", "a") as f:
        f.write(message1)
###################################################################################################
def analysis(file):
    QM_FMOPhore_analysis(file)
def FMOPhore_QM_analysis(pdb_files):
    QM_FMOPhore_analysis = analysis
    with Pool(processes=args.cpus) as pool:
        pool.map(QM_FMOPhore_analysis, pdb_files)
###################################################################################################    
def library_analysis(pdb_files):
    if args.same_target:
        pdb_processor = LibAnalysis(pdb_file=pdb_files, Binding_Energy=args.Binding_Energy)
        pdb_processor.Lib_analyze()
    else:
        pass
######################################## Message ###################################################
message3 = '''############################################################
# FMOPhore V.0.1 Copyright "©", Peter E.G.F. Ibrahim. 2024 #
############################################################\n'''
with open("FMOPhore.log", "a") as f:
    f.write(message3)
###################################################################################################
#  Timer 
###################################################################################################
def main_Processor():
    ###########################
    message_PDBProcessor = "PDBProcessor run..."
    message_PDBProcessor_fin = "Finished PDBProcessor successfully."
    with open("FMOPhore.log", "a") as log_file:
        log_file.write(message_PDBProcessor+"\n")
    preparePDBs_all(pdb_files)
    if args.align is True:
        message_align = "Alignment run..."
        message_align_fin = "Finished alignment successfully."
        with open("FMOPhore.log", "a") as log_file:
            log_file.write(message_align+"\n")
        if pdb_files:
            ref_file = pdb_files[0]
            pdb_files_notref = pdb_files[1:]
            align_library(ref_file, pdb_files_notref)
        with open("FMOPhore.log", "a") as log_file:
            log_file.write(message_align_fin+"\n")
    else:
        pass
    sep_LIGsPDBs_all(pdb_files)
    cutoff_PDBs_all(pdb_files)
    if args.same_target:
        lib_pdbs(pdb_files)
    else:
        pass
    splitPDBs_all(pdb_files)
    with open("FMOPhore.log", "a") as log_file:
        log_file.write(message_PDBProcessor_fin+"\n")
#################################################
def main_FMOPhore():
    ###########################
    message_FMOPhore = "FMOPhore run..."
    message_FMOPhore_fin = "Finished FMOPhore successfully."
    with open("FMOPhore.log", "a") as log_file:
        log_file.write(message_FMOPhore+"\n")
    FMOPhore_all(pdb_files)
#################################################
def main_QM():
    message_QM = '''QM processing...'''
    message_QM_fin = "Finished QM successfully."
    with open("FMOPhore.log", "a") as log_file:
        log_file.write(message_QM + "\n")
    QM(pdb_files, args.qm_calculation)
    with open("FMOPhore.log", "a") as log_file:
        log_file.write(message_QM_fin + "\n")
#################################################
def main_FMOPhore_analysis():
    message_QM_FMOPhoreAnalysis = "QM-FMOPhore Analysis processing..."
    message_QM_FMOPhoreAnalysis_fin = "Finished QM-FMOPhore Analysis successfully."
    with open("FMOPhore.log", "a") as log_file:
        log_file.write(message_QM_FMOPhoreAnalysis+"\n")
    FMOPhore_QM_analysis(pdb_files)
    with open("FMOPhore.log", "a") as log_file:
        log_file.write(message_QM_FMOPhoreAnalysis_fin+"\n")
#################################################
def main_library_analysis():
    message_QM_FMOPhoreAnalysis = "Library Analysis processing..."
    message_QM_FMOPhoreAnalysis_fin = "Finished Library Analysis successfully."
    with open("FMOPhore.log", "a") as log_file:
        log_file.write(message_QM_FMOPhoreAnalysis+"\n")
    library_analysis(pdb_files)
    with open("FMOPhore.log", "a") as log_file:
        log_file.write(message_QM_FMOPhoreAnalysis_fin+"\n")
#################################################
@timeout_decorator.timeout(args.timer * 86400)
def main():
    if args.qm_calculation:
        main_Processor()
        main_FMOPhore()
        print('''Please install GAMESS software: 
            https://www.msg.chem.iastate.edu/gamess/download.html
            If you have it installed and running:
            Recommended a GPU cluster equipped with at least 2 GPUs and 20 CPU.
            Contact the author Peter E.G.F. Ibrahim: 2448959@dundee.ac.uk, for further details on running FMOPhore to its full potential.''')

        # main_QM()
    else:
        pass
    time.sleep(10)    
    if args.qm_calculation or args.FMOPhore_analysis:
        main_FMOPhore_analysis()
        main_library_analysis()
    else:
        pass
    time.sleep(10)
##################################################################################################
##################################################################################################
if __name__ == '__main__':
    """
    FMOPhore V.0.1 - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    main()
