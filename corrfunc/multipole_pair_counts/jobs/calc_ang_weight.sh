#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=calc-ang-weight
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --output=../output/%x-%j.out

module load python/3.6

source ~/env/bin/activate

python ../code/calc_ang_weight.py $PAR_INPUT $FIB_INPUT $NORMS $OUTPUT

deactivate
