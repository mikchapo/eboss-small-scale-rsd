#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=dr-wp
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $DR_OUTPUT| ../code/05-comp_DR_eBOSS_proj_llist.exe
