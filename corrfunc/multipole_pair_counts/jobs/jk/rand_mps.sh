#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=rand-mps-jk
#SBATCH --output=../output/%x-%j.out

echo $RAND_CAT $RR_OUTPUT $NREG| ../code/13-comp_RR_eBOSS_llist_jk.exe
