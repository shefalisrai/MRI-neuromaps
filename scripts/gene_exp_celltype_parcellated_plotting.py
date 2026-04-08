#########VIZ########
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np


out_dir = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/'
measures = ['T1w_contrast', 'sulc', 'thickness', 'area', 'RND-gm', 'RNI-gm', 'RND-wm', 'RNI-wm']
contrasts = ['EDR_HC', 'EDBP_HC', 'EDR_EDBP']
celltypes_list = ['Astro', 'Endo', 'Micro', 'Neuro-Ex', 'Neuro-In', 'Oligo', 'OPC']


results_df = pd.read_csv(os.path.join(out_dir, 'celltype_correlation_results.csv'))

# Build r and p matrices: rows = measure|contrast, cols = celltypes
r_cols = [f'r_{c}' for c in celltypes_list]
p_cols = [f'p_{c}' for c in celltypes_list]

results_df['measure_contrast'] = results_df['measure'] + ' | ' + results_df['contrast']
results_df = results_df.set_index('measure_contrast')

pivot_r = results_df[r_cols].copy()
pivot_r.columns = celltypes_list

pivot_p = results_df[p_cols].copy()
pivot_p.columns = celltypes_list

sig_mask = pivot_p >= 0.05

#Plot
fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)

for ax, contrast in zip(axes, contrasts):
    rows = [f'{m} | {contrast}' for m in measures]
    r_sub = pivot_r.loc[rows].copy()
    p_sub = pivot_p.loc[rows].copy()
    mask_sub = sig_mask.loc[rows].copy()
    
    r_sub.index = measures
    p_sub.index = measures
    mask_sub.index = measures

    sns.heatmap(r_sub,
                cmap='RdBu_r',
                center=0,
                vmin=-0.6, vmax=0.6,  
                annot=r_sub.round(2),
                fmt='',
                linewidths=0.5,
                mask=mask_sub,
                ax=ax,
                cbar=False,           
                annot_kws={'size': 9})

    ax.set_title(contrast, fontsize=12, fontweight='bold')
    ax.set_ylabel('Measure', fontsize=10)
    ax.set_xlabel('')
    ax.tick_params(axis='y', labelsize=9)
    if ax == axes[-1]:
        ax.tick_params(axis='x', rotation=30, labelsize=10)

# Add single colorbar on the right
import matplotlib as mpl
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])  
norm = mpl.colors.Normalize(vmin=-0.6, vmax=0.6)
sm = mpl.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
sm.set_array([])
fig.colorbar(sm, cax=cbar_ax, label='Pearson r')

plt.suptitle('Cell type expression correlations\n(only p<0.05 shown)',
             fontsize=14, fontweight='bold')
plt.savefig(os.path.join(out_dir, 'celltype_heatmap.png'), dpi=150, bbox_inches='tight')
plt.close()
print('Saved celltype_heatmap.png')