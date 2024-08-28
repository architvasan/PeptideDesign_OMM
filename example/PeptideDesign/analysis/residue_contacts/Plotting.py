import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('contacts_1.csv')

df_idr = df[(df['resid']>119) & (df['resid']<192)]
df_tot = pd.DataFrame({'Metastate 2': list(df_idr['contactprob'])},
                        index=df_idr['resid'])

df = pd.read_csv(f'contacts_4.csv')
df_idr = df[(df['resid']>119) & (df['resid']<192)]
df_tot[f'Metastate 4'] = list(df_idr['contactprob'])

df_tot.to_csv('all_contact_probs_nmnat2_disord.csv', index=False)
print(df)

# Plot the DataFrame
ax = df_tot.plot(kind='bar', stacked=True, color=['skyblue','magenta'], figsize=(30, 10))

# Increase width of each bar
for container in ax.containers:
    for bar in container:
        bar.set_width(0.9)  # Adjust the width as needed


# Get the current x-tick positions
xticks = ax.get_xticks()

# Filter x-tick positions to show every 5th one
step = 5
filtered_xticks = xticks[::step]

# Set new x-ticks and format labels
ax.set_xticks(filtered_xticks)

# Set new x-ticks and format labels
ax.set_xticks(filtered_xticks)
#ax.set_xticklabels([f'{int(df_tot["resid"][int(tick)])}' for tick in filtered_xticks])#[f'{int(tick):,}' for tick in filtered_xticks])  # Format labels with commas

# Add title and labels
plt.xlabel('Resids', fontsize=14)
plt.ylabel('Contact Probability', fontsize=14)
plt.xticks(rotation=45)

# Show the plot
plt.savefig(f'contacts_nmnat2_disord_metas.png', bbox_inches='tight', dpi=300)
plt.close()
