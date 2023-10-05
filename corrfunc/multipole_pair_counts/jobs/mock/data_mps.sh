#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=data-mps-mock
#SBATCH --output=../output/%x-%j.out

echo $DATA_CAT $DD_OUTPUT| ../code/17-comp_DD_eBOSS_llist_mock.exe
