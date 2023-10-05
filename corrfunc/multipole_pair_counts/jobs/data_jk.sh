#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=data-jk-200-ngc-mock-0001
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --output=../output/%x-%j.out

echo ../data/EZmock_eBOSS_LRG_NGC_v7_0001_Nz-matched_regions_200.dat ../output/DD_EZmock_eBOSS_LRG_NGC_v7_0001_mps_jk_200.dat 201| ../code/11-comp_DD_eBOSS_llist_jk.exe
