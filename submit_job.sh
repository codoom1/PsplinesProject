#!/usr/bin/bash
#SBATCH --partition=cpu      # Partition (queue) name
#SBATCH --ntasks=1                   # Number of task
#SBATCH --cpus-per-task=6             # number of CPU cores per task
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --mem=40gb                    # Job memory request
#SBATCH --qos=long                      ## allows for more time
#SBATCH --time=4-00:00:00              # Time limit hrs:min:sec
#SBATCH --output=logs/parallel_r2_fast_%j.log    # Standard output and error log
#SBATCH --mail-type=BEGIN
#SBATCH --job-name=r2fastsim

module load r-rocker-ml-verse/4.4.0+apptainer

cd /PATH/TO/PROJECT
Rscript runL2rateCase.R