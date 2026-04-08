# -*- coding: utf-8 -*-
"""
Cell type gene expression correlations with spatial null models
Uses abagen cell type expression mapped to surface giftis not as csv
Then uses neuromaps for spin-test corrected correlations between cell type expression and 3 group contrasts for 8 measures.

SR, April 2026
"""

import os
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from neuromaps import datasets, resampling, stats, images, nulls

os.environ['PATH'] += ':/Applications/workbench/bin_macosx64'

gifti_dir    = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/3Groupcontrasts_funcgiftis/'
out_dir      = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/'
annot_lh     = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/parcellationDir/lh.Schaefer2018_100Parcels_7Networks_order.annot'
annot_rh     = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/parcellationDir/rh.Schaefer2018_100Parcels_7Networks_order.annot'
celltype_exp_path = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/schaefer100_celltype_expression.csv'
gifti_out_dir = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/celltype_giftis/'
os.makedirs(out_dir, exist_ok=True)
os.makedirs(gifti_out_dir, exist_ok=True)

measures  = ['T1w_contrast', 'sulc', 'thickness', 'area', 'RND-gm', 'RNI-gm', 'RND-wm', 'RNI-wm']
contrasts = ['EDR_HC', 'EDBP_HC', 'EDR_EDBP']
celltypes = ['Astro', 'Endo', 'Micro', 'Neuro-Ex', 'Neuro-In', 'Oligo', 'OPC']

#Annotations 
print('Loading annot files...')
lh_labels, _, _ = nib.freesurfer.read_annot(annot_lh)  # (163842,)
rh_labels, _, _ = nib.freesurfer.read_annot(annot_rh)  # (163842,)

print('Loading cell type expression...')
celltype_df = pd.read_csv(celltype_exp_path, index_col=0)
celltype_df = celltype_df.sort_index()  
print(f'Cell type expression shape: {celltype_df.shape}')

print('\nMapping cell type expression onto surface...')
celltype_annots = {}  

for ctype in celltypes:
    if ctype not in celltype_df.columns:
        print(f'  Skipping {ctype} - not found in expression matrix')
        continue

    vals = celltype_df[ctype].values  # (100,)

    # Map onto LH surface (parcels 1-50)
    lh_surface = np.zeros(163842, dtype=np.float32)
    for pid in range(1, 51):
        lh_surface[lh_labels == pid] = vals[pid - 1]

    # Map onto RH surface (parcels 51-100)
    rh_surface = np.zeros(163842, dtype=np.float32)
    for pid in range(1, 51):
        rh_surface[rh_labels == pid] = vals[pid + 49]

    # Save as gifti
    lh_path = os.path.join(gifti_out_dir, f'lh.{ctype}_celltype.func.gii')
    rh_path = os.path.join(gifti_out_dir, f'rh.{ctype}_celltype.func.gii')

    lh_img = nib.gifti.GiftiImage()
    lh_img.add_gifti_data_array(nib.gifti.GiftiDataArray(lh_surface))
    nib.save(lh_img, lh_path)

    rh_img = nib.gifti.GiftiImage()
    rh_img.add_gifti_data_array(nib.gifti.GiftiDataArray(rh_surface))
    nib.save(rh_img, rh_path)

    celltype_annots[ctype] = (lh_path, rh_path)
    print(f'  Saved {ctype} gifti')


print('\nResampling cell type giftis to fsaverage 10k...')

# Use genepc1 as target (fsaverage 10k)
genepc = datasets.fetch_annotation(desc='genepc1', return_single=True)

celltype_annots_10k = {}

for ctype, (lh_path, rh_path) in celltype_annots.items():
    annot_10k, _ = resampling.resample_images(  # swap _ and annot_10k
        src=(lh_path, rh_path),
        trg=genepc,
        src_space='fsaverage',
        trg_space='fsaverage',
        resampling='transform_to_trg')
    celltype_annots_10k[ctype] = annot_10k 
    print(f'  Resampled {ctype}')

print('\nRunning neuromaps correlations with spin tests...')

all_results = []

for measure in measures:
    for contrast in contrasts:
        lh_file = os.path.join(gifti_dir, f'lh.postpsubs_vertexwise_{contrast}_{measure}_sm256_baseline_3Group.zstat.func.gii')
        rh_file = os.path.join(gifti_dir, f'rh.postpsubs_vertexwise_{contrast}_{measure}_sm256_baseline_3Group.zstat.func.gii')

        if not os.path.exists(lh_file) or not os.path.exists(rh_file):
            print(f'[Skipping] {measure} - {contrast}')
            continue

        print(f'\nProcessing: {measure} - {contrast}')

        # Resample to fsaverage 10k
        my_map_10k, _ = resampling.resample_images(
            src=(lh_file, rh_file),
            trg=genepc,
            src_space='fsaverage',
            trg_space='fsaverage',
            resampling='transform_to_trg')

        # spins
        data_10k = images.load_data(my_map_10k)
        rot = nulls.alexander_bloch(data_10k,
                                     atlas='fsaverage', density='10k',
                                     n_perm=1000, seed=1234)

        row = {'measure': measure, 'contrast': contrast}

        #Correlate
        for ctype, annot_10k in celltype_annots_10k.items():
            r_val, p_val = stats.compare_images(my_map_10k, annot_10k, nulls=rot)
            row[f'r_{ctype}'] = r_val
            row[f'p_{ctype}'] = p_val
            print(f'  {ctype}: r={r_val:.3f}, p={p_val:.4f}')

        all_results.append(row)


results_df = pd.DataFrame(all_results)
results_df.to_csv(os.path.join(out_dir, 'celltype_correlation_spintests.csv'), index=False)
print(f'\nResults saved to {out_dir}celltype_correlation_spintests.csv')

####Plotting#####
r_cols = [f'r_{c}' for c in celltypes]
p_cols = [f'p_{c}' for c in celltypes]

results_df['measure_contrast'] = results_df['measure'] + ' | ' + results_df['contrast']
results_df = results_df.set_index('measure_contrast')

pivot_r = results_df[r_cols].copy()
pivot_r.columns = celltypes
pivot_p = results_df[p_cols].copy()
pivot_p.columns = celltypes
sig_mask = pivot_p >= 0.05

fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)

for ax, contrast in zip(axes, contrasts):
    rows = [f'{m} | {contrast}' for m in measures]
    r_sub = pivot_r.loc[rows].copy()
    p_sub = pivot_p.loc[rows].copy()
    mask_sub = sig_mask.loc[rows].copy()

    r_sub.index   = measures
    p_sub.index   = measures
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

#Onecolorbar at side
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
norm = mpl.colors.Normalize(vmin=-0.6, vmax=0.6)
sm = mpl.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
sm.set_array([])
fig.colorbar(sm, cax=cbar_ax, label='Pearson r')

plt.suptitle('Cell type expression correlations — spin test corrected\n(only p_spin < 0.05 shown)',
             fontsize=14, fontweight='bold')
plt.savefig(os.path.join(out_dir, 'celltype_heatmap_spintests.png'), dpi=150, bbox_inches='tight')
plt.close()
print('Saved celltype_heatmap_spintests.png')