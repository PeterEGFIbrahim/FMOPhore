# FMOPhore V.0.1 Hotspot prediction and classification

---

# <img width="900" alt="image" align="center" src="https://github.com/user-attachments/assets/4a3fbc8c-fd40-4b96-a621-dd14d669c0a3">

---

# FMOPhore v0.1 - Hotspot Prediction and Classification

FMOPhore is a Python package designed for hotspot identification and classification. It facilitates the preparation and processing of protein and ligand structures, as well as advanced quantum mechanics (QM) calculations. The tool supports various features, including binding energy calculations, and analysis MD trajectory processing.

---

## System Requirements

- **GAMESS Software**  
  Required for QM-FMO calculations.
  - GAMESS software (https://www.msg.chem.iastate.edu/gamess/download.html) to run QM-FMO calculations.
	Quantum Mechanics (QM) Calculations: To run QM calculations, ensure GAMESS is installed and contact the author for additional guidance.
  - [Download GAMESS](https://www.msg.chem.iastate.edu/gamess/download.html).  
  - A GPU cluster with at least 2 GPUs and 20 CPUs is recommended for QM calculations.  
  - Contact the author for additional guidance on running FMOPhore with GAMESS.
  - Please install GAMESS software: 
        If you have it installed and running:
        Recommended a GPU cluster equipped with at least 2 GPUs and 20 CPU.
        Contact the author Peter E.G.F. Ibrahim: pibrahim001@dundee.ac.uk - 2448959@dundee.ac.uk, for further details on running FMOPhore to its full potential.
    
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
- **Trajectory Support**: Perform MD using SuMD and ACEMD, for compatibility with FMOPhore analysis.
- **Analysis Tools**: Perform post-calculation analysis.

---
### Set up Conda Environment

To ensure all dependencies are correctly installed and avoid conflicts (Except: GAMESS, SuMD and ACEMD, mentioned above), follow these steps:

1. **Install Conda**: If Conda is not already installed, follow the installation guide [here](https://docs.conda.io/en/latest/miniconda.html).

2. **Create the Environment**:
   ```bash
   conda env create -f FMOPhore_env.yml
   conda activate FMOPhore_env
   ```
   
## Installation

### Prerequisites

- Python >= 3.10
- GAMESS software (for QM calculations).

### Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/PeterEGFIbrahim/FMOPhore.git

2. Navigate to the directory:
   ```bash
   cd FMOPhore
   
4. Install the package:
   ```bash
   pip install .
   
## Usage

### Mandatory Parameters:
   ```bash
    -dir,       --directory             : Process all PDB files in a directory
    -com,       --Prot_complex          : Process a single complex PDB file .pdb
    -pdb,       --PDB_ID                : PDB ID or .txt file containing PDB IDs.
                                          Process a single complex PDB directly from Protein Data Bank > https://www.rcsb.org/ 
    -prot       --protein_pdb_file      : Path to protein PDB file
    -ligs       --ligand_files          : Path to single ligand PDB file or directory of ligand PDB files

    -qm         --qm_calculation        : MP2 or DFTB (GAMESS software required/ https://www.msg.chem.iastate.edu/gamess/download.html)
    
    -d          --distance_cutoff       : Distance cutoff for selecting residues. Use "no_cutoff" to select the whole protein.
    
    -t          --timer                 : Timer in days e.g.: 1 day ==> -t 1  
    -c          --cpus                  : Please specify the number of cpus (cores; number of jobs to parallelize)
   ```
### Optional Parameters:
   ```bash
    -DA-FMO,    --trajectory            : To do DA-FMO, read trajectory "traj.dcd" 
                                        : Dy-FMOPhore requires SuMD (https://github.com/molecularmodelingsection/SuMD) 
                                        :                      ACEMD (https://software.acellera.com/acemd/tutorial.html)
    -cof,       --cofactor              : If you have co-factor in the system, provide the name in this form, e.g.: LYS-600'
    -BE,        --Binding_Energy        : Calculate the Binding Energy (ΔE) = (E_Complex) − (E_Protein) − (E_Ligand)
    -lib,       --same_target           : If same target and different ligands
    -align,     --align                 : If needed to align structures
    -analysis   --FMOPhore_analysis     : To run analysis only if calculations has been completed.
    -p          --personalized          : To run personalized
   ```


### Example Command
   ```
fmophore -dir /path/to/pdb/files -d 5 -qm DFTB -t 1 -c 10 -lib
   ```
### Help
To view the full list of options and their usage:
   ```
fmophore --help
   ```
### Outputs

After running FMOPhore, you can expect the following outputs:

- **FP_score.csv**: Contains the FP-score and interaction data.
- **FP_score_vs_Hotspots.png**: A plot visualizing FP scores vs. hotspots.
- **Ph4_heatmap.png**: A 2D heatmap of hotspots.

### Analysis Script

- Use `FP_score.py` to analyze and generate results.
- In Examples folder >> FP_score examples.

#### Sample Outputs:

- **merged_files.csv**: Contains Percentage of Interaction (%) and Lowest Energy per residue.
- **selected_data.csv**: Provides further processed results for detailed analysis.

### Requirements

The following Python libraries are required (automatically installed with the package):

- "numpy",
- "tqdm",
- "timeout-decorator",
- "argparse",
- "pandas",           
- "seaborn",
- "matplotlib",
- "rdkit"
  
### Developer Information

Author: Peter E.G.F. Ibrahim  
Email: pibrahim001@dundee.ac.uk, 2448959@dundee.ac.uk, peteregfi@gmail.com  
GitHub: [PeterEGFIbrahim](https://github.com/PeterEGFIbrahim)  


### Licensing
- This project is licensed under the GPL-3.0 License.
- See the LICENSE file for details.
    
### Copyright
© 2024 Peter E.G.F. Ibrahim. All rights reserved.
