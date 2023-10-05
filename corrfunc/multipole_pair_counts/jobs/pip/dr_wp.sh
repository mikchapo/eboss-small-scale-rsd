#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=dr-wp-pip-ang
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $BW_INPUT $DR_AUW_INPUT $DR_OUTPUT| ../code/10-comp_DR_eBOSS_proj_llist_pip.exe
