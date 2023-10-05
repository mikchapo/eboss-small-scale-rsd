#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=dr-mps-pip-ang
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $BW_INPUT $DR_AUW_INPUT $DR_OUTPUT| ../code/08-comp_DR_eBOSS_llist_pip.exe
