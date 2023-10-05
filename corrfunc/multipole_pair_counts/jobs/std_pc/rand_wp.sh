#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=rand-wp
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --output=../output/%x-%j.out

echo $RAND_INPUT $RR_OUTPUT| ../code/06-comp_RR_eBOSS_proj_llist.exe
