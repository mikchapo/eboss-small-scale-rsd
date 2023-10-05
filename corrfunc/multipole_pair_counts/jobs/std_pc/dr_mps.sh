#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=dr-mps
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $DR_OUTPUT| ../code/02-comp_DR_eBOSS_llist.exe
