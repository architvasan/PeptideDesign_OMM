# Peptide design project

## General steps
 1. Identify literature derived binding proteins to your protein of interest using PathwayCommons
 2. Generate multimer structure using ESMFold (can also use AlphaFold2 but saw better results with ESMFold)
 3. Create solvated solvated + ionized set up using ambertools (0.15 M NaCl)
 4. Simulate system using OpenMM for multiple replica
 5. Analyze coordinates of interest during simulation (e.g. ANCA to analyze IDRs)
 6. Construct MSM to identify metastable protein conformations
 7. Obtain conformations of protein of interest
 8. For each protein of interest, follow RFDiffusion protocol to generate novel binding peptides.

