#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=norm-count-RR-pip
#SBATCH --output=../output/%x-%j.out

echo $RAND_INPUT $NORM_RR_OUTPUT| ../code/27-norm_count_RR_pip.exe
