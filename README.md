# FMOPhore V.0.1 Hotspot prediction and classification

# <img width="900" alt="image" align="center" src="https://github.com/user-attachments/assets/4a3fbc8c-fd40-4b96-a621-dd14d669c0a3">

##   System requirments
- GAMESS software (https://www.msg.chem.iastate.edu/gamess/download.html) to run QM-FMO calculations.
	Quantum Mechanics (QM) Calculations: If you plan to run QM calculations, ensure GAMESS is installed and contact the author for additional guidance.
	System Requirements: For QM calculations, a GPU cluster equipped with at least 2 GPUs and 20 CPUs is recommended.
- SuMD (https://github.com/molecularmodelingsection/SuMD) and ACEMD (https://software.acellera.com/acemd/tutorial.html) to run Dy-FMOPhore analysis.

# FMOPhore v0.1 - Hotspot Prediction and Classification

FMOPhore is a Python package designed for hotspot identification and classification. It facilitates the preparation and processing of protein and ligand structures, as well as advanced quantum mechanics (QM) calculations. The tool supports various features, including binding energy calculations, trajectory processing, and analysis of previously completed computations.

---

## System Requirements

- **GAMESS Software**  
  Required for QM-FMO calculations.  
  - [Download GAMESS](https://www.msg.chem.iastate.edu/gamess/download.html).  
  - A GPU cluster with at least 2 GPUs and 20 CPUs is recommended for QM calculations.  
  - Contact the author for additional guidance on running FMOPhore with GAMESS.

- **SuMD and ACEMD**  
  Required for Dy-FMOPhore analysis.  
  - [SuMD](https://github.com/molecularmodelingsection/SuMD)  
  - [ACEMD](https://software.acellera.com/acemd/tutorial.html)

---

## Features

- **Hotspot Identification and Classification**: Identify and classify hotspots using advanced QM methods.
- **PDB Preparation**: Process and prepare PDB files for calculations.
- **QM Calculations**: Run MP2 or DFTB calculations using GAMESS software.
- **Binding Energy Calculations**: Calculate the binding energy (ΔE) for protein-ligand complexes.
- **Trajectory Support**: Analyze trajectories with tools like SuMD and ACEMD.
- **Analysis Tools**: Perform post-calculation analysis.

---

## Installation

### Prerequisites

- Python >= 3.6
- GAMESS software (optional, for QM calculations).

### Set up Conda Environment

To ensure all dependencies are correctly installed and avoid conflicts, follow these steps:

1. **Install Conda**: If Conda is not already installed, follow the installation guide [here](https://docs.conda.io/en/latest/miniconda.html).

2. **Create the Environment**:
   ```bash
   conda env create -f FMOPhore_env.yml
   conda activate FMOPhore_env


### Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/PeterEGFIbrahim/FMOPhore.git

Navigate to the directory:
   ```bash
cd FMOPhore

Install the package:
   ```bash
pip install .


### Mandatory Parameters:
   ```bash
-dir,       --directory             : Process all PDB files in a directory.
-com,       --Prot_complex          : Process a single complex PDB file.
-PDB,       --PDB_ID                : Specify a PDB ID or a file containing multiple PDB IDs.
-prot,      --protein_pdb_file      : Path to the protein PDB file.
-ligs,      --ligand_files          : Path to a single ligand PDB file or a directory of ligand PDB files.

-PDBProcessor, --PDBProcessor       : Prepare the PDBs only.
    -d         --distance_cutoff    : Distance cutoff for selecting residues (or "no_cutoff" for the whole protein).

-FMOPhore, --FMOPhore               : Run FMOPhore.
    -qm       --qm_calculation      : Specify "MP2" or "DFTB" (requires GAMESS software).

-t,         --timer                 : Timer in days (e.g., `-t 1` for 1 day).
-c,         --cpus                  : Number of CPUs to use for parallelization.

### Optional Parameters:
   ```bash
-DA-FMO,    --trajectory            : Analyze trajectories (requires SuMD and ACEMD).
-cof,       --cofactor              : Specify co-factor (e.g., `LYS-600`).
-BE,        --Binding_Energy        : Calculate binding energy (ΔE).
-lib,       --same_target           : Specify if analyzing the same target with different ligands.
-align,     --align                 : Align structures if needed.
-analysis,  --FMOPhore_analysis     : Perform analysis of completed calculations.

Example Command

fmophore -dir /path/to/pdb/files -d 5 -qm DFTB -t 1 -c 20 -PDBProcessor -FMOPhore -lib -align

Help
To view the full list of options and their usage:

fmophore --help

Outputs
After running FMOPhore, you can expect the following outputs:

FP_score.csv: Contains the FP-score and interaction data.
FP_score_vs_Hotspots.png: A plot visualizing FP scores vs. hotspots.
Ph4_heatmap.png: A 2D heatmap of hotspots.

Analysis Script
Use FP_score.py to analyze and generate results.
Sample Outputs:

merged_files.csv and selected_data.csv (Percentage of Interaction (%), Lowest Energy per residue).

Requirements
The following Python libraries are required (automatically installed with the package):

numpy
tqdm
timeout-decorator
argparse

Developer Information
Author: Peter E.G.F. Ibrahim
Email: 2448959@dundee.ac.uk, peteregfi@gmail.com
GitHub: PeterEGFIbrahim

Licensing
This project is licensed under the GPL-3.0 License. See the LICENSE file for details.
    
Copyright
© 2024 Peter E.G.F. Ibrahim. All rights reserved.
