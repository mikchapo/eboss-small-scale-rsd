export DATA_INPUT_NGC="../data/eBOSS_LRG_NGC_pip_v7_2.dat"
export BW_INPUT_NGC="../data/eBOSS_LRG_NGC_pip_v7_2_bw.dat"
export DATA_INPUT_SGC="../data/eBOSS_LRG_SGC_pip_v7_2.dat"
export BW_INPUT_SGC="../data/eBOSS_LRG_SGC_pip_v7_2_bw.dat"
export NORM_DD_CROSS_OUTPUT="../output/DD_norm_cross_eBOSS_LRG_v7_2_pip.dat"
#sbatch pip/norm_count_DD_cross.sh
#for cap in "NGC" "SGC"
for cap in "SGC"
  do
  export DATA_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2.dat"
  export BW_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2_bw.dat"
  export DD_AUW_INPUT="../output/DD_ang_weights_eBOSS_LRG_${cap}_v7_2.dat"
  export DR_AUW_INPUT="../output/DR_ang_weights_eBOSS_LRG_${cap}_v7_2.dat"
  export RAND_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2.rand"
  export NORM_DD_OUTPUT="../output/DD_norm_eBOSS_LRG_${cap}_v7_2_pip.dat"
  export NORM_DR_OUTPUT="../output/DR_norm_eBOSS_LRG_${cap}_v7_2_pip.dat"
  export NORM_RR_OUTPUT="../output/RR_norm_eBOSS_LRG_${cap}_v7_2_pip.dat"
#  sbatch pip/norm_count_DD.sh
#  sbatch pip/norm_count_DR.sh
#  sbatch pip/norm_count_RR.sh
#  for cf in "mps" "wp"
  for cf in "wp"
    do
    export DD_OUTPUT="../output/DD_eBOSS_LRG_${cap}_v7_2_${cf}_pip_ang_test.dat"
    export DR_OUTPUT="../output/DR_eBOSS_LRG_${cap}_v7_2_${cf}_pip_ang.dat"
    export RR_OUTPUT="../output/RR_eBOSS_LRG_${cap}_v7_2_${cf}_pip.dat"
    sbatch pip/data_${cf}.sh
#    sbatch pip/dr_${cf}.sh
#    sbatch pip/rand_${cf}.sh
    done
  done


