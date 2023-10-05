#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=analyze-chains-INPUT_NAME
#SBATCH --output=../output/job_logs/analyze_chains/%x-%j.out

module load scipy-stack
source ~/P2/fit/emEnv/bin/activate

python ../code/analyze_chains_P2.py $INPUT_FILE

deactivate
