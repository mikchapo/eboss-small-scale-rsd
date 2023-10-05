for cap in "NGC" "SGC" 
  do
  export DATA_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2_new.dat.fits"
  export RAND_INPUT="../data/eBOSS_LRG_${cap}_v7_2_pip_new.ran.fits"
  export DATA_OUTPUT="../output/eBOSS_LRG_${cap}_pip_v7_2.dat"
  export BW_OUTPUT="../output/eBOSS_LRG_${cap}_pip_v7_2_bw.dat"
  export RAND_OUTPUT="../output/eBOSS_LRG_${cap}_pip_v7_2.rand"
#  sbatch assemble_data.sh
  sbatch assemble_rand.sh
  done
