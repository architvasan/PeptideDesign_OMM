source ~/.bashrc
conda activate /lambda_stor/homes/avasan/miniforge3/envs/deeptime

frame_file_list=("metafrmes_0_0.dat" "metafrmes_4_0.dat" "metafrmes_0_1.dat" "metafrmes_0_2.dat" "metafrmes_3_2.dat" "metafrmes_0_3.dat" "metafrmes_1_3.dat" "metafrmes_5_3.dat" "metafrmes_0_4.dat" "metafrmes_1_4.dat" "metafrmes_2_4.dat")
trajn_list=("0" "0" "1" "2" "2" "3" "3" "3" "4" "4" "4")


# Loop over both arrays at the same time
    echo "${array1[$i]} is ${array2[$i]}"
#done


for ((i=0; i<${#frame_file_list[@]}; i++)); do

    echo "Processing $file"
    python ../../src/analyze/msm/metastable_trajs.py -p ../Protein_OnlySims/protein_rep0.pdb -t ../Protein_OnlySims/prot_rep${trajn_list[$i]}.dcd -F out/meta/${frame_file_list[$i]} -O out/meta -ot meta${frame_file_list[$i]}.dcd
    # Perform actions on $file
done

