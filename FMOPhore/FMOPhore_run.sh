#!/usr/bin/env bash
#$ -pe smp 10
#$ -jc long
#$ -cwd
#$ -p -100
#qsub -l gpu=0 -mods m_mem_free 16G

conda activate FMOPhore_env
 
# Holo-FMOPhore
fmophore -dir ./ -d 5 -qm DFTB -t 1 -c 10 -lib

# Apo-scan-FMOPhore
# fmophore -prot prot/*.pdb -ligs ligands/ -d 5 -qm DFTB -t 1 -c 10 -lib 

conda deactivate
conda deactivate