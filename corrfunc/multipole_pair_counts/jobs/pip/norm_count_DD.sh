#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=norm-count-DD-pip
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $BW_INPUT $NORM_DD_OUTPUT| ../code/23-norm_count_DD_pip.exe
