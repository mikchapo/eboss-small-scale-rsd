#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=norm-count-DR-pip
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $BW_INPUT $NORM_DR_OUTPUT| ../code/26-norm_count_DR_pip.exe
