for cap in "NGC" "SGC"
  do
  for pc in "DD" "DR"
    do
    export PAR_INPUT="../output/${pc}_eBOSS_LRG_${cap}_v7_2_ang_par.dat"
    export FIB_INPUT="../output/${pc}_eBOSS_LRG_${cap}_v7_2_ang_fib_pip.dat"
    export NORMS="../output/${pc}_norm_eBOSS_LRG_${cap}_v7_2_ang.dat"
    export OUTPUT="../output/${pc}_ang_weights_eBOSS_LRG_${cap}_v7_2.dat"
    sbatch calc_ang_weight.sh
    done
  done


