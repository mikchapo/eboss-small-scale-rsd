#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=rand-jk-100-ngc
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --output=../output/%x-%j.out

echo ../archive/jk_100/data/eBOSS_LRG_clustering_NGC_v7_hinv_regions_100.rand ../archive/jk_100/output/RR_eBOSS_LRG_clustering_NGC_v7_rs_llist_ab_jk_100_sparse.dat 101| ../code/13-comp_RR_eBOSS_llist_jk.exe
