#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=data-mps
#SBATCH --output=../output/%x-%j.out

echo $DATA_INPUT $DD_OUTPUT| ../code/01-comp_DD_eBOSS_llist.exe
