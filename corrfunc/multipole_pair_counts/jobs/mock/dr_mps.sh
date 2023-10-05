#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=dr-mps-mock
#SBATCH --output=../output/%x-%j.out

echo $DATA_CAT $RAND_CAT $DR_OUTPUT| ../code/18-comp_DR_eBOSS_llist_mock.exe
