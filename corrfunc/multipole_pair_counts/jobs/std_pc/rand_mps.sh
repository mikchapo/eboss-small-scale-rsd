#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=rand-mps
#SBATCH --output=../output/%x-%j.out

echo $RAND_INPUT $RR_OUTPUT| ../code/03-comp_RR_eBOSS_llist.exe
