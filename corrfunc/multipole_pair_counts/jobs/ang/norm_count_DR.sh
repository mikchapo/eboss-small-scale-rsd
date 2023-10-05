#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=norm-count-DR-ang
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $BW_INPUT $NORM_DR_OUTPUT| ../code/29-norm_count_DR_ang.exe
