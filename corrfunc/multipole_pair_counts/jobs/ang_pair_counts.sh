for cap in "NGC" "SGC"
  do
  export DATA_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2.dat"
  export BW_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2_bw.dat"
  export RAND_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2.rand"
  export NORM_DD_OUTPUT="../output/DD_norm_eBOSS_LRG_${cap}_v7_2_ang.dat"
  export NORM_DR_OUTPUT="../output/DR_norm_eBOSS_LRG_${cap}_v7_2_ang.dat"
  sbatch ang/norm_count_DD.sh
  sbatch ang/norm_count_DR.sh
  export DD_PAR_OUTPUT="../output/DD_eBOSS_LRG_${cap}_v7_2_ang_par.dat"
  export DR_PAR_OUTPUT="../output/DR_eBOSS_LRG_${cap}_v7_2_ang_par.dat"
  export DD_FIB_OUTPUT="../output/DD_eBOSS_LRG_${cap}_v7_2_ang_fib_pip.dat"
  export DR_FIB_OUTPUT="../output/DR_eBOSS_LRG_${cap}_v7_2_ang_fib_pip.dat"
#  sbatch ang/data.sh
#  sbatch ang/dr.sh
  done


