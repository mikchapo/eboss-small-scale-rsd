#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=dr-ang
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $BW_INPUT $DR_PAR_OUTPUT $DR_FIB_OUTPUT| ../code/25-comp_DR_eBOSS_ang.exe
