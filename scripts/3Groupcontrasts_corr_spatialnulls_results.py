# -*- coding: utf-8 -*-
"""
Neuromaps spatial null model correlations
EDC clinic vertexwise data vs neuromaps annotations

Using 3 group contrasts (EDR-HC, EDBP-HC, EDR-EDBP) and 8 measures (T1w_contrast, sulc, thickness, area, RND-gm, RNI-gm, RND-wm, RNI-wm)
**n_perm = 100 (for quick testing right now, but should be 1000+ for final analysis)**

March 26, 2026
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
out_dir = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/'
os.makedirs(out_dir, exist_ok=True)

measures = ['T1w_contrast', 'sulc', 'thickness', 'area', 'RND-gm', 'RNI-gm', 'RND-wm', 'RNI-wm']
contrasts = ['EDR_HC', 'EDBP_HC', 'EDR_EDBP']

#Annotations-fsaverage 10k
genepc = datasets.fetch_annotation(desc='genepc1', return_single=True) #gene PC1

#Annotations-fsLR 164k
subjvar = datasets.fetch_annotation(desc='intersubjvar', return_single=True) #intersubject variability
cbf = datasets.fetch_annotation(source='raichle', desc='cbf', return_single=True) #cerebral blood flow
cbv = datasets.fetch_annotation(source='raichle', desc='cbv', return_single=True) #cerebral blood volume
cmrglc_raichle = datasets.fetch_annotation(source='raichle', desc='cmrglc', return_single=True) #glucose metabolism
cmro2 = datasets.fetch_annotation(source='raichle', desc='cmr02', return_single=True) #oxygen metabolism

#Specific DevExp annotation (right hemi only)
devexp = datasets.fetch_annotation(source='hill2010', desc='devexp', return_single=True)
rh_sphere = '/Users/shefalirai/neuromaps-data/atlases/fsLR/tpl-fsLR_den-164k_hemi-R_sphere.surf.gii'

#Annotations-fsLR 32k
saaxis = datasets.fetch_annotation(source='sydnor2021', desc='SAaxis', return_single=True) #sensorimotor-association axis
evoexp_xu = datasets.fetch_annotation(source='xu2020', desc='evoexp', return_single=True) #evolutionary expansion

#Grouped by the fs/fsLR spaces 
fsaverage_annots = [
 ('genepc1', genepc),
]
fslr164k_annots = [
 ('subjvar', subjvar),
 ('cbf', cbf),
 ('cbv', cbv),
 ('cmrglc_raichle', cmrglc_raichle),
 ('cmro2', cmro2)
]
fslr32k_annots = [
 ('SAaxis', saaxis),
 ('evoexp_xu', evoexp_xu)
]
 
#############Output Results#############
results = []
for measure in measures:
    for contrast in contrasts:
 
        lh = os.path.join(gifti_dir, f'lh.postpsubs_vertexwise_{contrast}_{measure}_sm256_baseline_3Group.zstat.func.gii')
        rh = os.path.join(gifti_dir, f'rh.postpsubs_vertexwise_{contrast}_{measure}_sm256_baseline_3Group.zstat.func.gii')
 
        if not os.path.exists(lh) or not os.path.exists(rh):
            print(f'[Skipping] Missing files for {measure} - {contrast}')
            continue
 
        my_map = (lh, rh)
        print(f'\nRunning: {measure} - {contrast}')
 
        row = {'measure': measure, 'contrast': contrast}
 
        #fsaverage 10k
        my_map_r1, _ = resampling.resample_images(
            src=my_map, trg=genepc,
            src_space='fsaverage', trg_space='fsaverage',
            resampling='transform_to_trg')
 
        data_r1 = images.load_data(my_map_r1)
        rot_r1  = nulls.alexander_bloch(data_r1, atlas='fsaverage', density='10k',
                                        n_perm=1000, seed=1234)
 
        for name, annot in fsaverage_annots:
            _, annot_r = resampling.resample_images(
                src=my_map, trg=annot,
                src_space='fsaverage', trg_space='fsaverage',
                resampling='transform_to_trg')
            r_val, p_val = stats.compare_images(my_map_r1, annot_r, nulls=rot_r1)
            row[f'r_{name}'] = r_val
            row[f'p_{name}'] = p_val
            print(f'  {name}: r = {r_val:.3f}, p = {p_val:.4f}')
 
        #fsLR 164k
        my_map_r2, _ = resampling.resample_images(
            src=my_map, trg=subjvar,
            src_space='fsaverage', trg_space='fsLR',
            resampling='transform_to_trg')
 
        data_r2 = images.load_data(my_map_r2)
        rot_r2  = nulls.alexander_bloch(data_r2, atlas='fsLR', density='164k',
                                        n_perm=1000, seed=1234)
 
        for name, annot in fslr164k_annots:
            _, annot_r = resampling.resample_images(
                src=my_map, trg=annot,
                src_space='fsaverage', trg_space='fsLR',
                resampling='transform_to_trg')
            r_val, p_val = stats.compare_images(my_map_r2, annot_r, nulls=rot_r2)
            row[f'r_{name}'] = r_val
            row[f'p_{name}'] = p_val
            print(f'  {name}: r = {r_val:.3f}, p = {p_val:.4f}')
 
        #Separate resampling for Devexp-RH only
        my_map_rh, _ = resampling.resample_images(
            src=rh, trg=devexp,
            src_space='fsaverage', trg_space='fsLR',
            resampling='transform_to_trg',
            hemi='R')

        data_rh = images.load_data(my_map_rh)
        rot_rh  = nulls.alexander_bloch(data_rh, atlas='fsLR', density='164k',
                                        n_perm=1000, seed=1234,
                                        surfaces=(rh_sphere,))

        _, annot_r = resampling.resample_images(
            src=rh, trg=devexp,
            src_space='fsaverage', trg_space='fsLR',
            resampling='transform_to_trg',
            hemi='R')

        r_val, p_val = stats.compare_images(my_map_rh, annot_r, nulls=rot_rh)
        row['r_devexp'] = r_val
        row['p_devexp'] = p_val
        print(f'  devexp: r = {r_val:.3f}, p = {p_val:.4f}')

        #fsLR 32k
        my_map_r3, _ = resampling.resample_images(
            src=my_map, trg=saaxis,
            src_space='fsaverage', trg_space='fsLR',
            resampling='transform_to_trg')
 
        data_r3 = images.load_data(my_map_r3)
        rot_r3  = nulls.alexander_bloch(data_r3, atlas='fsLR', density='32k',
                                        n_perm=1000, seed=1234)
 
        for name, annot in fslr32k_annots:
            _, annot_r = resampling.resample_images(
                src=my_map, trg=annot,
                src_space='fsaverage', trg_space='fsLR',
                resampling='transform_to_trg')
            r_val, p_val = stats.compare_images(my_map_r3, annot_r, nulls=rot_r3)
            row[f'r_{name}'] = r_val
            row[f'p_{name}'] = p_val
            print(f'  {name}: r = {r_val:.3f}, p = {p_val:.4f}')
 
        results.append(row)
 
#############SaveCSV#############
df = pd.DataFrame(results)
df.to_csv(os.path.join(out_dir, 'neuromaps_results.csv'), index=False)
print(f'\nResults saved to {out_dir}neuromaps_results.csv')
print(df.to_string())
 
