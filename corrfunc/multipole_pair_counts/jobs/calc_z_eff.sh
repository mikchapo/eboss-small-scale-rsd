#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-z-eff
#SBATCH --output=../output/%x-%j.out

echo  "../data/eBOSS_LRG_SGC_pip_v7_2.dat" "../data/eBOSS_LRG_SGC_pip_v7_2_bw.dat" ../output/z_eff_eBOSS_LRG_SGC_v7_2_pip.dat| ../code/32-calc_z_eff_pip.exe
