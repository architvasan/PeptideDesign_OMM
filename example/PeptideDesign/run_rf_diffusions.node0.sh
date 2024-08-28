#mkdir outputs_1
../src/peptide_design/run_rfdiffusion.sh inputpdbs/meta0_cleaned.pdb outputs_0/peptides A122,A162,A165,A166,A167,A169,A170,A172 30 307 20 100 log_meta0 0 /eagle/datascience/avasan/Simulations/NMNAT-2/Simulations_NMNAT-2_FBXO45/src/peptide_design &

#mkdir outputs_1
mkdir outputs_1_b1
../src/peptide_design/run_rfdiffusion.sh inputpdbs/meta1_cleaned.pdb outputs_1_b1/peptides A121,A122,A123 30 307 20 100 log_meta1_b1 1 /eagle/datascience/avasan/Simulations/NMNAT-2/Simulations_NMNAT-2_FBXO45/src/peptide_design &

mkdir outputs_1_b2
../src/peptide_design/run_rfdiffusion.sh inputpdbs/meta1_cleaned.pdb outputs_1_b2/peptides A162,A167,A169,A170,A171,A172 30 307 20 100 log_meta1_b2 2 /eagle/datascience/avasan/Simulations/NMNAT-2/Simulations_NMNAT-2_FBXO45/src/peptide_design &

#mkdir outputs_2
mkdir outputs_2_b1
../src/peptide_design/run_rfdiffusion.sh inputpdbs/meta2_cleaned.pdb outputs_2/peptides A121,A122,A123 30 307 20 100 log_meta2_b1 3 /eagle/datascience/avasan/Simulations/NMNAT-2/Simulations_NMNAT-2_FBXO45/src/peptide_design 

#../src/peptide_design/run_rfdiffusion.sh inputpdbs/meta2_cleaned.pdb outputs_2/peptides A121,A122,A123,A162,A165,A167,A169,A170,A171,A172 30 307 20 100 log_meta2 3 /eagle/datascience/avasan/Simulations/NMNAT-2/Simulations_NMNAT-2_FBXO45/src/peptide_design &

#mkdir outputs_4
#../src/peptide_design/run_rfdiffusion.sh inputpdbs/meta4_cleaned.pdb outputs_4/peptides A164,A165,A166,A167,A169,A170 30 307 20 100 log_meta4 3 /eagle/datascience/avasan/Simulations/NMNAT-2/Simulations_NMNAT-2_FBXO45/src/peptide_design

