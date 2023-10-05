for cap in "NGC" "SGC"
  do
  export DATA_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2.dat"
  export RAND_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2.rand"
  export NORM_OUTPUT="../output/norm_eBOSS_LRG_${cap}_v7_2_cp.dat"
  sbatch std_pc/norm_count.sh
#for cf in "mps" "wp"
for cf in "wp"
    do
    export DD_OUTPUT="../output/DD_eBOSS_LRG_${cap}_v7_2_${cf}_cp.dat"
    export DR_OUTPUT="../output/DR_eBOSS_LRG_${cap}_v7_2_${cf}_cp.dat"
    export RR_OUTPUT="../output/RR_eBOSS_LRG_${cap}_v7_2_${cf}_cp.dat"
#    sbatch std_pc/data_${cf}.sh
#    sbatch std_pc/dr_${cf}.sh
#    sbatch std_pc/rand_${cf}.sh
    done
  done


