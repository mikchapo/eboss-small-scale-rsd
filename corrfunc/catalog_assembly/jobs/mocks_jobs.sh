for mock_id in "0001" "0002" "0003" "0004" "0005" "0006" "0007" "0008" "0009" "0010" 
  do
  export MOCK_ID=$mock_id
  sbatch assemble_mock.sh
  done
