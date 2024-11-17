# FMOPhore

#####################################################################

•	Full automation software for fragment hit to lead growth and optimization: 

System preparation step:
	o	Preparation, fragmentation. 
	o	Automated Multichain pdbs processing.
	o	Automated alignment of complexes for complete library of molecules or snapshots during MD for DA-QM-FMO processing.

QM-DA-FMO/ QM-FMO calculation steps:
	o	Calculating ligand-protein binding interaction energies.
	o	Calculating ligand energy, protein energy and protein-ligand complex energy (full system).
	o	Binding Energy = (E complex) – (E protein) – (E ligand).
	o	Preforming Dynamical-averaging-FMO (DA-FMO) on a full Molecular dynamics trajectory of ligand-protein complex.
	o	QM-FMO calculations (DFTB-MP2) on protein, ligand and complex separate.

Binding site and hotspots indentificaiton steps:
	o	Descriptors generation. 
	o	Binding site and hotspot identification.

Pharmacophores and growth vectors steps:
	o	Holo-pharmacophore generation.
	o	Growth vectors prediction (Holo and Apo growth vectors). 

Analysis and plots generation step:
	o	Giving full details about the ligands 3D coordinates in the binding site, with the distances and angles of binding interaction, bond types (HB, hydrophobic, ionic, Pi-Pi stacking …etc) and interaction energies, binding energy, and deltaG, plots PEIDA and heatmaps.

Note: Make sure the system is an output from Schrodniger software preparation step. 

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
