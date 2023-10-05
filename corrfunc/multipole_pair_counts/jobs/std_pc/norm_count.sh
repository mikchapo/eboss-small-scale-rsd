#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=norm-count
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $NORM_OUTPUT| ../code/00-weight_count_eBOSS.exe
