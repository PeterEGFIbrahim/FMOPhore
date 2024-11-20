# FMOPhore

<img width="459" alt="image" src="https://github.com/user-attachments/assets/4a3fbc8c-fd40-4b96-a621-dd14d669c0a3">


#####################################################################

# Conda Environment installation
	conda activate
	conda env create -f FMOPhore_env.yml
	conda activate FMOPhore_env

# Default to run on the gpu-s directly

	python FraGrow_run.py -dir ./ -d 5 -lib same -qm DFTB &

# Helper
	python FraGrow_run.py -h

:||: FraGrow 1.0 :||:

optional arguments:
  -h, --help            show this help message and exit

Mandatory parameters:

    -dir,       --directory             : Process all PDB files in a directory
    -com,       --Prot_complex          : Process a single complex PDB file .pdb
    -PDB,       --PDB_ID                : PDB ID or .txt file containing PDB IDs.
                                          Process a single complex PDB directly from Protein Data Bank > https://www.rcsb.org/ 
    -prot       --protein_pdb_file      : Path to protein PDB file
    -ligs       --ligand_files          : Path to single ligand PDB file or directory of ligand PDB files

    -d          --distance_cutoff       : Distance cutoff for selecting residues. Use "no_cutoff" to select the whole protein.
    -qm         --qm_calculation        : MP2 or DFTB
    -t          --timer                 : Timer in days e.g.: 1 day ==> -t 1  
:
  -These are required parameters THESE ARE REQUIRED PARAMETERS, -! THESE ARE REQUIRED PARAMETERS

Optinal parameters:

    -DA-FMO,    --trajectory            : To do DA-FMO, read trajectory "traj.dcd"
    -cof,       --cofactor              : If you have co-factor in the system, provide the name in this form, e.g.: LYS-600'
    -lib,       --same_target           : If same target and different ligands
:
  -These are optinal parameters THESE ARE OPTINAL PARAMETERS, -... THESE ARE OPTINAL PARAMETERS

Peter E.G.F. Ibrahim.


# To run on the cluster using the queuing system 
# Open the FraGrow.sh file and specify the command-line.

	qsub FraGrow.sh 

#####################################################################

Copyright "©" Developer: Peter E. G. F. Ibrahim.
