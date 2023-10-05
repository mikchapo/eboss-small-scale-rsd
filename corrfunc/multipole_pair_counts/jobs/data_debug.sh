#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=data-debug-sgc-1c-1o-2g
#SBATCH --output=../output/%x-%j.out

echo ../data/eBOSS_LRG_clustering_SGC_v7.dat ../output/DD_debug.dat| ../code/05-comp_DD_eBOSS_llist_debug.exe
