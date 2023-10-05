#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=16G
#SBATCH --account=rr-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=sham-hod-check
#SBATCH --output=../output/job_logs/%x-%j.out

module load scipy-stack

cd ../code/
python sham_hod_check.py
