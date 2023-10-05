#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=data-wp-pip-ang
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $BW_INPUT $DD_AUW_INPUT $DD_OUTPUT| ../code/09-comp_DD_eBOSS_proj_llist_pip.exe
