import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df_pae = pd.DataFrame({})
for i in range(0, 6):
    
    pae_dict = {}
    if i<3:
        pae_dict[f'Metastable {i+1}'] = np.loadtxt(f'pae_int_{i}.dat')
    elif i==3:
        continue
    elif i>3:
        pae_dict[f'Metastable {i}'] = np.loadtxt(f'pae_int_{i}.dat')

    df_pae_it = pd.DataFrame(pae_dict)
    df_pae = pd.concat([df_pae, df_pae_it], ignore_index=True, axis=1)
    print(df_pae)
#print(pae_dict)
df_pae.plot(kind='box')
plt.xlabel('Metastable state')
plt.ylabel('PAE interaction')
plt.savefig(f'pae_ints.png', dpi=300, bbox_inches='tight')
plt.close()
