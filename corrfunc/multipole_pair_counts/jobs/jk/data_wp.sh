#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=data-wp-jk
#SBATCH --output=../output/%x-%j.out

echo $DATA_CAT $DD_OUTPUT $NREG| ../code/14-comp_DD_eBOSS_proj_llist_jk.exe
