#!/bin/bash
# Check if the correct number of arguments is provided
if [ "$#" -lt 5 ]; then
  echo "Usage: $0 path_mpnn path_af2 input_silent mpnn_out af2_out"
  exit 1
fi

path_mpnn=$1
path_af2=$2
input_silent=$3
mpnn_out=$4
af2_out=$5

module use /soft/modulefiles
module load conda
conda activate /lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design

echo "${size_ref}"

CUDA_VISIBLE_DEVICES=0,1,2,3 ${path_mpnn}/dl_interface_design.py -silent $input_silent -outsilent $mpnn_out

##################################### AlphaFold2 Complex Prediction ##########################################
##############################################################################################################

conda deactivate
conda activate /lus/eagle/projects/datascience/avasan/envs/af2_binder_design
${path_af2}/predict.py -silent $mpnn_out -outsilent $af2_out

