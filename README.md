# Peptide design project

## General steps
 1. Identify literature derived binding proteins to your protein of interest using PathwayCommons
    - https://apps.pathwaycommons.org/
 2. Generate multimer structure using ESMFold (can also use AlphaFold2 but saw better results with ESMFold)
    
    `python src/fold_ai/esm_pred.py -h`
 3. Create solvated solvated + ionized set up using ambertools (0.15 M NaCl)
    
    `python src/create/create_apo.py -h`
 4. Simulate system using OpenMM for multiple replica
    
    `python run_openmm_mpi.py -h`
 5. Analyze coordinates of interest during simulation (e.g. ANCA to analyze IDRs)
    
    `python src/analyze/heiranca_dihedrals.py -h`
 6. Construct MSM to identify metastable protein conformations
    
    ```python
    python src/analyze/msm/preprocess.py -h
    python src/analyze/msm/msm.py -h
    python src/analyze/msm/metastable_trajs.py -h

    ```

 7. ID residue hotspots for each conformation
    `python src/analyze/contacts_A_B.py -h`
 8. For each protein of interest conformation, follow RFDiffusion protocol to generate novel binding peptides.
    - Scripts in `src/peptide_design`
    - Follow instructions in https://github.com/RosettaCommons/RFdiffusion.git for installation + running

