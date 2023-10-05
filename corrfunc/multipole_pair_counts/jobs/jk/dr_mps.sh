#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=dr-mps-jk
#SBATCH --output=../output/%x-%j.out

echo $DATA_CAT $RAND_CAT $DR_OUTPUT $NREG| ../code/12-comp_DR_eBOSS_llist_jk.exe
