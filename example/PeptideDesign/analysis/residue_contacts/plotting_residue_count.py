import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#df_tot = pd.DataFrame({'Metastable state': [], 'Positives': [], 'Negatives': [], 'Polars': [], 'NonPolars': []})
for i in range(6): 
    if i==3:
        continue
    
    if i!=2 and i!=1:
        df = pd.read_csv(f'analysis/residue_dists/residue_dist_{i}.csv', usecols=['Positives', 'Negatives', 'Polars', 'NonPolars'])
        df = pd.concat([df, pd.read_csv(f'analysis/residue_dists/residue_dist_{i}_small.csv',usecols=['Positives', 'Negatives', 'Polars', 'NonPolars'])]) 

    if i==1 or i==2:
        df = pd.read_csv(f'analysis/residue_dists/residue_dist_{i}_b1.csv',usecols=['Positives', 'Negatives', 'Polars', 'NonPolars'])
        df = pd.concat([df, pd.read_csv(f'analysis/residue_dists/residue_dist_{i}_b2.csv',usecols=['Positives', 'Negatives', 'Polars', 'NonPolars'])])    
        df = pd.concat([df, pd.read_csv(f'analysis/residue_dists/residue_dist_{i}_b1_small.csv',usecols=['Positives', 'Negatives', 'Polars', 'NonPolars'])])
        df = pd.concat([df, pd.read_csv(f'analysis/residue_dists/residue_dist_{i}_b2_small.csv',usecols=['Positives', 'Negatives', 'Polars', 'NonPolars'])])

    df_norm = df.div(df.sum(axis=1), axis=0)
    df_norm['Positives'].plot.kde(color='blue')#kind='hist', bins=10)
    df_norm['Negatives'].plot.kde(color='red')#(kind='hist', bins=10)
    df_norm['Polars'].plot.kde(color='green')#(kind='hist', bins=10)
    df_norm['NonPolars'].plot.kde(color='gray')#(kind='hist', bins=10)
    plt.savefig(f'distribution_residues_{i}.png', dpi=300, bbox_inches='tight')
    plt.close()
    #if i<3:
    #    df_norm['Metastable state'] = [i+1 for it in range(len(df_norm))]
    #else:
    #    df_norm['Metastable state'] = [i for it in range(len(df_norm))]

    #df_tot = pd.concat([df_tot, df_norm])
    #df_norm.to_csv(f'analysis/residue_dists/probs_res_{i}.csv', index=False)
    
