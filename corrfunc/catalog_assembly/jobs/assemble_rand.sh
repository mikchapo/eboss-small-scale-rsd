#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=assemble-rand-pip
#SBATCH --output=../output/%x-%j.out

module load python/3.6

source ~/env/bin/activate

python ../code/assemble_rand.py $RAND_INPUT $RAND_OUTPUT

deactivate
