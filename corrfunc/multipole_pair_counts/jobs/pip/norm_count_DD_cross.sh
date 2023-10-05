#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=norm-count-DD-pip-cross
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT_NGC $DATA_INPUT_SGC $BW_INPUT_NGC $BW_INPUT_SGC $NORM_DD_CROSS_OUTPUT| ../code/30-norm_count_DD_pip_cross.exe
