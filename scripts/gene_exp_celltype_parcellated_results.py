# -*- coding: utf-8 -*-
"""
Cell type gene expression correlations using manual pearson correlations (no spatial nulls!)
Here we have abagen cell type expression not mapped to surface giftis, but parcellated and saved as a csv file
Then we correlate the parcellated cell type expression with the parcellated group contrast maps (also parcellated in a separate script) using simple pearson correlations.

SR, April 2026
"""

import os
import numpy as np
import pandas as pd
import nibabel as nib
from scipy import stats
from statsmodels.stats.multitest import multipletests
from neuromaps import resampling, datasets, images

gifti_dir = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/3Groupcontrasts_funcgiftis/'
out_dir = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/'
annot_lh = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/parcellationDir/lh.Schaefer2018_100Parcels_7Networks_order.annot'
annot_rh = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/parcellationDir/rh.Schaefer2018_100Parcels_7Networks_order.annot'
gene_exp_path = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/schaefer100_gene_expression.csv'
celltype_path = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/celltypes_PSP.csv'
os.makedirs(out_dir, exist_ok=True)

measures = ['T1w_contrast', 'sulc', 'thickness', 'area', 'RND-gm', 'RNI-gm', 'RND-wm', 'RNI-wm']
contrasts = ['EDR_HC', 'EDBP_HC', 'EDR_EDBP']

# Load annot labels
print('Loading annot files...')
lh_labels, _, _ = nib.freesurfer.read_annot(annot_lh)
rh_labels, _, _ = nib.freesurfer.read_annot(annot_rh)

# Load full gene expression (100 parcels x 15633 genes)
print('Loading gene expression...')
gene_exp = pd.read_csv(gene_exp_path, index_col=0)
gene_exp = gene_exp.sort_index()  # ensure 1-100 order

# Load cell type gene lists and compute cell type expression
print('Computing cell type expression...')
celldata = pd.read_csv(celltype_path)
print(f'Cell types: {celldata["class"].unique()}')
celltypes = celldata['class'].unique()

celltype_exp = {}
for ctype in celltypes:
    geneid = celldata.loc[celldata['class'] == ctype, 'gene']
    ind = gene_exp.columns.isin(geneid)
    n_found = ind.sum()
    print(f'  {ctype}: {n_found}/{len(geneid)} genes found in expression matrix')
    if n_found > 0:
        genedata = gene_exp.iloc[:, ind].values  # (100, n_genes)
        celltype_exp[ctype] = np.mean(genedata, axis=1)  # (100,) mean across genes

celltype_df = pd.DataFrame(celltype_exp)  # (100 parcels x n_celltypes)
celltype_df.index = gene_exp.index
print(f'\nCell type expression shape: {celltype_df.shape}')
print(celltype_df.head())

# Save cell type expression for reference
celltype_df.to_csv(os.path.join(out_dir, 'schaefer100_celltype_expression.csv'))

# Load fsaverage 164k annotation for resampling
sert = datasets.fetch_annotation(source='beliveau2017', desc='dasb',
                                  space='fsaverage', den='164k',
                                  return_single=True)

parcel_ids_lh = np.unique(lh_labels[lh_labels > 0])  # 1-50
parcel_ids_rh = np.unique(rh_labels[rh_labels > 0])  # 1-50

all_results = []

for measure in measures:
    for contrast in contrasts:
        lh_file = os.path.join(gifti_dir, f'lh.postpsubs_vertexwise_{contrast}_{measure}_sm256_baseline_3Group.zstat.func.gii')
        rh_file = os.path.join(gifti_dir, f'rh.postpsubs_vertexwise_{contrast}_{measure}_sm256_baseline_3Group.zstat.func.gii')

        if not os.path.exists(lh_file) or not os.path.exists(rh_file):
            print(f'[Skipping] {measure} - {contrast}')
            continue

        print(f'\nProcessing: {measure} - {contrast}')

        # Resample from fsaverage 10k -> fsaverage 164k
        my_map_164k, _ = resampling.resample_images(
            src=(lh_file, rh_file),
            trg=sert,
            src_space='fsaverage',
            trg_space='fsaverage',
            resampling='transform_to_trg')

        all_data = images.load_data(my_map_164k)
        lh_data = all_data[:163842]
        rh_data = all_data[163842:]

        # Parcellate
        parcel_means = []
        for pid in parcel_ids_lh:
            mask = lh_labels == pid
            vals = lh_data[mask]
            parcel_means.append(np.mean(vals) if len(vals) > 0 else np.nan)
        for pid in parcel_ids_rh:
            mask = rh_labels == pid
            vals = rh_data[mask]
            parcel_means.append(np.mean(vals) if len(vals) > 0 else np.nan)

        my_map = np.array(parcel_means)  # (100,)
        valid_mask = ~np.isnan(my_map)

        # Correlate with each cell type
        row = {'measure': measure, 'contrast': contrast}
        for ctype in celltype_df.columns:
            ctype_vals = celltype_df[ctype].values
            r, p = stats.pearsonr(my_map[valid_mask], ctype_vals[valid_mask])
            row[f'r_{ctype}'] = r
            row[f'p_{ctype}'] = p
            print(f'  {ctype}: r={r:.3f}, p={p:.4f}')

        all_results.append(row)

# Save
results_df = pd.DataFrame(all_results)
results_df.to_csv(os.path.join(out_dir, 'celltype_correlation_results.csv'), index=False)
print(f'\nDone! Saved to {out_dir}celltype_correlation_results.csv')
print(results_df.to_string())
