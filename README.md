# FMOPhore

# <img width="900" alt="image" align="center" src="https://github.com/user-attachments/assets/4a3fbc8c-fd40-4b96-a621-dd14d669c0a3">


#####################################################################
<!-- System requirments: -->
- GAMESS software (https://www.msg.chem.iastate.edu/gamess/download.html)
- SuMD (https://github.com/molecularmodelingsection/SuMD)
- ACEMD (https://software.acellera.com/acemd/tutorial.html)
#####################################################################
<!-- Hardware recommended and tested: -->
- Recommended a GPU cluster equipped with at least 2 GPUs and 20 CPU.
- Tested on 2 NVIDIA GTX 1080 per node on 8 GPUs.
#####################################################################
Outputs from GAMESS software and PLIP, are combined in the FP-score equation.
 as shown in the sample attached. 
 merged_files.csv and selected_data.csv (shows the Percentage of Interaction (%), Lowest Energy per residue)
-Expected outputs:
 1. FP_score.csv 
 2. FP_score_vs_Hotspots.png (FP-score plot)
 3. Ph4_heatmap.png (2D-FMOPhore heatmaps)
#####################################################################
# Conda Environment installation
	conda activate
	conda env create -f FMOPhore_env.yml
	conda activate FMOPhore_env
# run
	python FP_score.py

Copyright "©" Developer: Peter E. G. F. Ibrahim.
