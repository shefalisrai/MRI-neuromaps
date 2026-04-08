# -*- coding: utf-8 -*-
"""
Neuromaps viz script
Loads results CSV and plots without rerunning the full analysis
March 26, 2026
SR
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

out_dir = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/'
df = pd.read_csv(out_dir + 'neuromaps_results.csv')
print(f'Loaded {len(df)} rows')
print(df.to_string())

contrast_order = ['EDR_HC', 'EDBP_HC', 'EDR_EDBP']
measure_order  = ['T1w_contrast', 'sulc', 'thickness', 'area', 'RND-gm', 'RNI-gm', 'RND-wm', 'RNI-wm']
all_annot_names = [c.replace('r_', '') for c in df.columns if c.startswith('r_')]

#############VISUALIZE#############
plot_rows = []
for _, row in df.iterrows():
    for name in all_annot_names:
        plot_rows.append({
            'annotation':  name,
            'measure':     row['measure'],
            'contrast':    row['contrast'],
            'r':           row[f'r_{name}'],
            'p':           row[f'p_{name}'],
            'abs_r':       abs(row[f'r_{name}']),
            'significant': row[f'p_{name}'] < 0.05,
        })
plot_df = pd.DataFrame(plot_rows)

#Individual heatmaps
for name in all_annot_names:
    sub     = plot_df[plot_df['annotation'] == name]
    r_pivot = sub.pivot(index='measure', columns='contrast', values='r').reindex(index=measure_order, columns=contrast_order)
    p_pivot = sub.pivot(index='measure', columns='contrast', values='p').reindex(index=measure_order, columns=contrast_order)

    annot_cells = r_pivot.copy().astype(object)
    for i in r_pivot.index:
        for j in r_pivot.columns:
            sig = '*' if p_pivot.loc[i, j] < 0.05 else ''
            annot_cells.loc[i, j] = f'{r_pivot.loc[i, j]:.2f}{sig}'

    fig, ax = plt.subplots(figsize=(7, 4))
    sns.heatmap(r_pivot, annot=annot_cells, fmt='', cmap='RdBu_r', center=0,
                vmin=-0.5, vmax=0.5, linewidths=0.5, ax=ax, annot_kws={'size': 8})
    ax.set_title(f'Neuromaps correlation: {name}\n(* = p < 0.05, spin test)')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='x', labelsize=8, rotation=30)
    ax.tick_params(axis='y', labelsize=8, rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'heatmap_{name}.png'), dpi=150)
    plt.close()
    print(f'Saved: heatmap_{name}.png')

#FULL UNIFIEDheatmap--per contrast, annotations x measures
fig, axes = plt.subplots(1, 3, figsize=(len(all_annot_names) * 1.8, len(measure_order) * 0.85 + 1.5),
                         sharey=True)

for ax_idx, (ax, contrast) in enumerate(zip(axes, contrast_order)):
    sub     = plot_df[plot_df['contrast'] == contrast]
    r_pivot = sub.pivot(index='measure', columns='annotation', values='r').reindex(
                  index=measure_order, columns=all_annot_names)
    p_pivot = sub.pivot(index='measure', columns='annotation', values='p').reindex(
                  index=measure_order, columns=all_annot_names)

    annot_cells = r_pivot.copy().astype(object)
    for i in r_pivot.index:
        for j in r_pivot.columns:
            sig = '*' if p_pivot.loc[i, j] < 0.05 else ''
            annot_cells.loc[i, j] = f'{r_pivot.loc[i, j]:.2f}{sig}'

    sns.heatmap(r_pivot, annot=annot_cells, fmt='', cmap='RdBu_r', center=0,
                vmin=-0.5, vmax=0.5, linewidths=0.5, ax=ax,
                cbar=ax_idx == 2,
                annot_kws={'size': 7})

    ax.set_title(contrast, fontsize=11, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='x', labelsize=10, rotation=45)
    ax.tick_params(axis='y', labelsize=10, rotation=0)

fig.suptitle('Spatial Correlations\n(* = p < 0.05, spin test)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'unified_heatmap_neuromaps.png'), dpi=500, bbox_inches='tight')
plt.close()
print('Saved: unified_heatmap_neuromaps.png')