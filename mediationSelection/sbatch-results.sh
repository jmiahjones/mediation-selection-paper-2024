#!/bin/bash
#SBATCH --partition=standard -J med-nopass
#SBATCH -c 24 --mem=250G
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_90
#SBATCH -t 5-00:00
#SBATCH --array=0-23
#SBATCH -d after:8502215_15

echo Array index: $SLURM_ARRAY_TASK_ID
hostname
date

module load r/3.6.3/b1

cores=24 ./run-sims-ml-sbatch-innerpar.sh
