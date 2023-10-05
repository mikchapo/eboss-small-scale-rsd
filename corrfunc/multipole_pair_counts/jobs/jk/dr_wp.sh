#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=dr-wp-jk
#SBATCH --output=../output/%x-%j.out

echo $DATA_CAT $RAND_CAT $DR_OUTPUT $NREG| ../code/15-comp_DR_eBOSS_proj_llist_jk.exe
