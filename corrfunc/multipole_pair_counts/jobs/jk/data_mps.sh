#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=data-mps-jk
#SBATCH --output=../output/%x-%j.out

echo $DATA_CAT $DD_OUTPUT $NREG| ../code/11-comp_DD_eBOSS_llist_jk.exe
