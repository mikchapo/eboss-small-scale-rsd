#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=rand-wp-jk
#SBATCH --mail-user=mj3chapm@uwaterloo.ca      
#SBATCH --mail-type=ALL
#SBATCH --output=../output/%x-%j.out

echo $RAND_CAT $RR_OUTPUT $NREG| ../code/16-comp_RR_eBOSS_proj_llist_jk.exe
