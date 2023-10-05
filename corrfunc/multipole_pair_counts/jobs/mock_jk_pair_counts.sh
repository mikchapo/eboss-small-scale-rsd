for mock_id in "0001" "0002" "0003" "0004" "0005" "0006" "0007" "0008" "0009" "0010" 
# for mock_id in "0001"
  do
  export MOCK_ID=$mock_id
  export NREG=201
  for cap in "NGC" "SGC"
    do
    export DATA_CAT="../data/EZmock_eBOSS_LRG_${cap}_v7_${mock_id}_regions_200.dat"
    export RAND_CAT="../data/EZmock_eBOSS_LRG_${cap}_v7_${mock_id}_regions_200.rand"
    export NORM_OUTPUT="../output/norms_EZmock_eBOSS_LRG_${cap}_v7_${mock_id}.dat"
    sbatch jk/weight_count.sh
    for cf in "mps" "wp"
      do
      export DD_OUTPUT="../output/DD_EZmock_eBOSS_LRG_${cap}_v7_${mock_id}_${cf}_jk_200.dat"
      export DR_OUTPUT="../output/DR_EZmock_eBOSS_LRG_${cap}_v7_${mock_id}_${cf}_jk_200.dat"
      export RR_OUTPUT="../output/RR_EZmock_eBOSS_LRG_${cap}_v7_${mock_id}_${cf}_jk_200.dat"
#      sbatch jk/data_${cf}.sh
#      sbatch jk/dr_${cf}.sh
#      sbatch jk/rand_${cf}.sh
      done
    done
  done
