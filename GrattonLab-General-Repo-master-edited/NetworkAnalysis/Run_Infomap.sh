#!/bin/bash
#SBATCH -A b1081
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -t 100:00:00
#SBATCH -p b1081
#SBATCH --mem=0G
#SBATCH --job-name="infomap"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ally.dworetsky@northwestern.edu
#SBATCH -o "%x.o%j"


## job commands; <matlabscript> is your MATLAB .m file, specified without the .$
module load matlab/r2020b
matlab -nosplash -nodesktop -r "addpath(genpath('/projects/b1081/Scripts/CIFTI_RELATED')); Run_Infomap_GrattonLab('"$1"'); quit"
