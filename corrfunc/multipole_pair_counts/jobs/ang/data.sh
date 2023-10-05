#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=data-ang
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $BW_INPUT $DD_PAR_OUTPUT $DD_FIB_OUTPUT| ../code/24-comp_DD_eBOSS_ang.exe
