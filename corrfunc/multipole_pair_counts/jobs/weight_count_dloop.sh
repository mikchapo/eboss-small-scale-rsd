#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=weight-count-dloop-uw
#SBATCH --output=../output/%x-%j.out

echo ../data/eBOSS_ELG_clustering_eboss22_v4_weighted.cat ../data/eBOSS_ELG_clustering_eboss22_v4_rand.cat ../output/eBOSS_ELG_clustering_eboss22_weight_counts_dloop_unweighted.txt unweighted| ../code/08-weight_count_dloop.exe
