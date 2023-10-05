#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=256M
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=assemble-mock
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --output=../output/%x-%j.out

module load python/3.6

source ~/env/bin/activate

python ../code/assemble_mock.py $MOCK_ID

deactivate
