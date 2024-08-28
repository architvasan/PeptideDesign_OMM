import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

for it in range(6):
    resid_data = pd.read_csv(f'out/meta/contacts_meta{it}/prob_nmnat2_disord.dat')
    # Create a bar plot
    plt.bar(resid_data['resid #'], resid_data[' contact prob'], color='blue')
    
    # Add title and labels
    plt.xlabel('Resids')
    plt.ylabel('Contact Probability')
    
    # Show the plot
    plt.savefig(f'images/contacts_nmnat2_disord_meta{it}.png', bbox_inches='tight', dpi=300)
    plt.close()
