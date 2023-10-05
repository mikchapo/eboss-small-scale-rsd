for input_name in cmass+eboss_baseline cmass+eboss_large_only eboss_ao_baseline eboss_ao_probes eboss_ao_scales_v3 eboss_cosmo_priors sham_ao
  do
  export INPUT_FILE="/home/mj3chapm/RSD/fit/input/tp_inputs/${input_name}.yaml"
  sed "s/INPUT_NAME/${input_name}/" analyze_chains_template.sh > analyze_chains_run.sh
  sbatch analyze_chains_run.sh
  sleep 1
  done
