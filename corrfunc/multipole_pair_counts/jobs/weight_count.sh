#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=weight-count-sgc
#SBATCH --output=../output/%x-%j.out

echo ../archive/mock_0001/data/EZmock_eBOSS_LRG_SGC_v7_0001_Nz-matched_mike.dat ../archive/mock_0001/data/EZmock_eBOSS_LRG_SGC_v7_0001_Nz-matched_mike.rand  ../output/norms_EZmock_eBOSS_LRG_SGC_v7_0001.txt| ../code/04-weight_count_eBOSS.exe
