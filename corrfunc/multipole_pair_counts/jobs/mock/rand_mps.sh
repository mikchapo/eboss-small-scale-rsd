#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=rand-mps-mock
#SBATCH --output=../output/%x-%j.out

echo $RAND_CAT $RR_OUTPUT| ../code/19-comp_RR_eBOSS_llist_mock.exe
