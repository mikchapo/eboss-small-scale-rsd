#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=dr-wp-mock
#SBATCH --output=../output/%x-%j.out

echo $DATA_CAT $RAND_CAT $DR_OUTPUT| ../code/21-comp_DR_eBOSS_proj_llist_mock.exe
