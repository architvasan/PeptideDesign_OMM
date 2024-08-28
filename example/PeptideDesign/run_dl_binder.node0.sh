#!/bin/bash
#mkdir outputs_1

path_mpnn=/lus/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/mpnn_fr
path_af2=/lus/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/af2_initial_guess

../src/peptide_design/run_dl_binder.sh $path_mpnn $path_af2 ./outputs_0/peptides.silent outputs_0/peptides_mpnn_out.silent outputs_0/peptides_af2_out.silent 0 logs/dl_0 &

../src/peptide_design/run_dl_binder.sh $path_mpnn $path_af2 ./outputs_1_b1/peptides.silent outputs_1_b1/peptides_mpnn_out.silent outputs_1_b1/peptides_af2_out.silent 1 logs/dl_1b1 &

../src/peptide_design/run_dl_binder.sh $path_mpnn $path_af2 ./outputs_1_b2/peptides.silent outputs_1_b2/peptides_mpnn_out.silent outputs_1_b2/peptides_af2_out.silent 2 logs/dl_1b2 &

../src/peptide_design/run_dl_binder.sh $path_mpnn $path_af2 ./outputs_2/peptides.silent outputs_2/peptides_mpnn_out.silent outputs_2/peptides_af2_out.silent 3 logs/dl_2 

