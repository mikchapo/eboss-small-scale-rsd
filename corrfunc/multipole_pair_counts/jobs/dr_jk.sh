#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=16G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=dr-jk-200-sgc-mock-0001
#SBATCH --output=../output/%x-%j.out

echo ../data/EZmock_eBOSS_LRG_SGC_v7_0001_Nz-matched_regions_200.dat ../data/EZmock_eBOSS_LRG_SGC_v7_0001_Nz-matched_regions_200.rand ../output/DR_EZmock_eBOSS_LRG_SGC_v7_0001_mps_jk_200.dat 201| ../code/12-comp_DR_eBOSS_llist_jk.exe
