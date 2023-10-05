for cap in "NGC" "SGC"
  do
  export DATA_CAT="../data/EZmock_eBOSS_LRG_${cap}_v7_0001_Nz-matched_mike.dat"
  export RAND_CAT="../data/EZmock_eBOSS_LRG_${cap}_v7_0001_Nz-matched_faizan.rand"
  for cf in "mps" "wp"
    do
    export DD_OUTPUT="../output/DD_EZmock_eBOSS_LRG_${cap}_v7_0001_${cf}_nc_zr_fr.dat"
    export DR_OUTPUT="../output/DR_EZmock_eBOSS_LRG_${cap}_v7_0001_${cf}_nc_zr_fr.dat"
    export RR_OUTPUT="../output/RR_EZmock_eBOSS_LRG_${cap}_v7_0001_${cf}_nc_zr_fr.dat"
    sbatch mock/data_${cf}.sh
    sbatch mock/dr_${cf}.sh
    sbatch mock/rand_${cf}.sh
    done
  done


