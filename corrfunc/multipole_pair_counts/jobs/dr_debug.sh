#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=dr-ngc-debug
#SBATCH --output=../output/%x-%j.out

echo ../data/eBOSS_LRG_clustering_NGC_v7.dat ../data/eBOSS_LRG_clustering_NGC_v7.rand ../output/DR_eBOSS_LRG_clustering_NGC_v7_rs_llist_debug.dat| ../code/06-comp_DR_eBOSS_llist_debug.exe
