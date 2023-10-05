#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=norm-count-DD-ang
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $BW_INPUT $NORM_DD_OUTPUT| ../code/28-norm_count_DD_ang.exe
