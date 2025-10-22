#!/bin/bash
#SBATCH -p compute
#SBATCH -n 25
#SBATCH --mem 100GB
#SBATCH --nodes=1
#SBATCH -t 0-24:00 # time (D-HH:MM)
#SBATCH -o slurm.loop.%j.out
#SBATCH -e slurm.loop.%j.err
#SBATCH --array=1-40%8  # Adjust 1-N to #lines in celltype_list.txt

###################
module load singularity

# read cell type list (tab-separated, no header; assumes one per line)
cell_types=($(tail -n +$SLURM_ARRAY_TASK_ID ./input/avannan_celltype_list.txt | head -1))

# escape slashes 
id_new=$(echo ${cell_types[0]} | sed 's_/_\\/_g')

# calc proximity
singularity exec ./containers/geospatial_latest.sif Rscript ./code/av_cellproximity_calculation_share.R "$id_new"

