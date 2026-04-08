# -*- coding: utf-8 -*-
"""
Neuromaps spatial null model correlations
EDC clinic vertexwise data vs neuromaps annotations

Using 3 group contrasts (EDR-HCxage, EDBP-HCxage, EDR-EDBPxage) and 8 measures
Comparing against xu2020 evoexp (bilateral) using groupxage contrasts

March 2026
SR
"""

import os
os.environ['PATH'] += ':/Applications/workbench/bin_macosx64'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from neuromaps import datasets, resampling, stats, images, nulls

gifti_dir = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/3Groupcontrasts_funcgiftis/'
out_dir   = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/'
os.makedirs(out_dir, exist_ok=True)

measures  = ['T1w_contrast', 'sulc', 'thickness', 'area', 'RND-gm', 'RNI-gm', 'RND-wm', 'RNI-wm']
contrasts = ['EDR_HC', 'EDBP_HC', 'EDR_EDBP'] #actually groupxage interactions, but named this way for simplicity

# xu2020 evoexp — bilateral (L and R)
evoexp_xu = datasets.fetch_annotation(source='xu2020', desc='evoexp')  # returns (L, R)

fslr32k_annots = [
    ('evoexp_xu', evoexp_xu),
]

###Output Results###
results = []

for measure in measures:
    for contrast in contrasts:

        lh = os.path.join(gifti_dir, f'lh.postpsubs_vertexwise_{contrast}_{measure}_sm256_baseline_3Group_interaction.zstat.func.gii')
        rh = os.path.join(gifti_dir, f'rh.postpsubs_vertexwise_{contrast}_{measure}_sm256_baseline_3Group_interaction.zstat.func.gii')

        if not os.path.exists(lh) or not os.path.exists(rh):
            print(f'[Skipping] Missing files for {measure} - {contrast}')
            continue

        print(f'\nRunning: {measure} - {contrast}')
        row = {'measure': measure, 'contrast': contrast}

        my_map = (lh, rh)

        # Resample to fsLR 32k (xu2020 space)
        my_map_r, _ = resampling.resample_images(
            src=my_map,
            trg=evoexp_xu,
            src_space='fsaverage',
            trg_space='fsLR',
            resampling='transform_to_trg'
        )

        # Spin nulls — bilateral fsLR 32k
        data_r = images.load_data(my_map_r)
        rot_r = nulls.alexander_bloch(
            data_r,
            atlas='fsLR',
            density='32k',
            n_perm=1000,
            seed=1234
        )

        for name, annot in fslr32k_annots:
            _, annot_r = resampling.resample_images(
                src=my_map,
                trg=annot,
                src_space='fsaverage',
                trg_space='fsLR',
                resampling='transform_to_trg'
            )
            r_val, p_val = stats.compare_images(my_map_r, annot_r, nulls=rot_r)
            row[f'r_{name}'] = r_val
            row[f'p_{name}'] = p_val
            print(f'  {name}: r = {r_val:.3f}, p = {p_val:.4f}')

        results.append(row)

###Save CSV###
df = pd.DataFrame(results)
df.to_csv(os.path.join(out_dir, 'neuromaps_evoexp_bilateral_results.csv'), index=False)
print(f'\nResults saved to {out_dir}neuromaps_evoexp_bilateral_results.csv')
print(df.to_string())

# ###Viz###
# all_annot_names = [name for name, _ in fslr32k_annots]
# contrast_order  = ['EDR_HC', 'EDBP_HC', 'EDR_EDBP']
# measure_order   = ['T1w_contrast', 'sulc', 'thickness', 'area', 'RND-gm', 'RNI-gm', 'RND-wm', 'RNI-wm']

# plot_rows = []
# for _, row in df.iterrows():
#     for name in all_annot_names:
#         plot_rows.append({
#             'annotation':  name,
#             'measure':     row['measure'],
#             'contrast':    row['contrast'],
#             'r':           row[f'r_{name}'],
#             'p':           row[f'p_{name}'],
#             'significant': row[f'p_{name}'] < 0.05,
#         })
# plot_df = pd.DataFrame(plot_rows)

# for name in all_annot_names:
#     sub     = plot_df[plot_df['annotation'] == name]
#     r_pivot = sub.pivot(index='measure', columns='contrast', values='r').reindex(index=measure_order, columns=contrast_order)
#     p_pivot = sub.pivot(index='measure', columns='contrast', values='p').reindex(index=measure_order, columns=contrast_order)

#     annot_cells = r_pivot.copy().astype(object)
#     for i in r_pivot.index:
#         for j in r_pivot.columns:
#             sig = '*' if p_pivot.loc[i, j] < 0.05 else ''
#             annot_cells.loc[i, j] = f'{r_pivot.loc[i, j]:.2f}{sig}'

#     fig, ax = plt.subplots(figsize=(7, 4))
#     sns.heatmap(r_pivot, annot=annot_cells, fmt='', cmap='RdBu_r', center=0,
#                 vmin=-0.5, vmax=0.5, linewidths=0.5, ax=ax, annot_kws={'size': 8})
#     ax.set_title(f'Neuromaps correlation: {name} (bilateral)\n(* = p < 0.05, spin test)')
#     ax.set_xlabel('')
#     ax.set_ylabel('')
#     ax.tick_params(axis='x', labelsize=8, rotation=30)
#     ax.tick_params(axis='y', labelsize=8, rotation=0)
#     plt.tight_layout()
#     plt.savefig(os.path.join(out_dir, f'heatmap_bilateral_{name}.png'), dpi=150)
#     plt.close()
#     print(f'Saved: heatmap_bilateral_{name}.png')