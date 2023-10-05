#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=data-mps-pip-ang
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $BW_INPUT $DD_AUW_INPUT $DD_OUTPUT| ../code/07-comp_DD_eBOSS_llist_pip.exe
