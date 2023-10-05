#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=assemble-data-pip
#SBATCH --output=../output/%x-%j.out

module load python/3.6

source ~/env/bin/activate

python ../code/assemble_data.py $DATA_INPUT $DATA_OUTPUT $BW_OUTPUT

deactivate
