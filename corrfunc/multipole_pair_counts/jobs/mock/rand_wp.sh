#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=rand-wp-mock
#SBATCH --output=../output/%x-%j.out

echo $RAND_CAT $RR_OUTPUT| ../code/22-comp_RR_eBOSS_proj_llist_mock.exe
