#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=weight-count-mock
#SBATCH --output=../output/%x-%j.out

echo $DATA_CAT $RAND_CAT $NORM_OUTPUT| ../code/04-weight_count_eBOSS.exe
