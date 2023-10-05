#!/bin/bash
#SBATCH --time=02:30:00
#SBATCH --cpus-per-task=32
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=rand-wp-pip
#SBATCH --output=../output/%x-%j.out

echo $RAND_INPUT $RR_OUTPUT| ../code/06-comp_RR_eBOSS_proj_llist.exe
