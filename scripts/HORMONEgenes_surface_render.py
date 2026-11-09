# -*- coding: utf-8 -*-
"""
Render hormone/ED candidate-gene AHBA expression on the cortical surface
(lateral + medial, magma), matching the cell-type gene-expression figure.

Robust two-stage approach: render each (gene, view) to its own small PNG with
nilearn (fast), then composite the grid with imshow (instant). Avoids packing
many 3D axes into one figure, which hangs.

Maps: AHBA (abagen) Schaefer-100 expression -> surface giftis in hormonegene_giftis/.
"""
import os
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from nilearn import datasets, plotting

base      = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps'
gifti_dir = os.path.join(base, 'hormonegene_giftis')
out_dir   = os.path.join(base, 'results')
tmp_dir   = os.path.join(out_dir, '_hormone_tiles')
os.makedirs(tmp_dir, exist_ok=True)

genes = ['ESR1', 'GPER1', 'AR', 'OXTR', 'LEPR', 'BDNF', 'FTO']
gene_labels = {'ESR1':'ESR1 (ERα)', 'GPER1':'GPER1', 'AR':'AR',
               'OXTR':'OXTR', 'LEPR':'LEPR', 'BDNF':'BDNF', 'FTO':'FTO'}

fs = datasets.fetch_surf_fsaverage('fsaverage')
infl_l = fs['infl_left']
cmap = 'magma'

# ---- stage 1: render each gene lateral + medial (LH) to small PNGs ----
for g in genes:
    lh = nib.load(os.path.join(gifti_dir, f'lh.{g}_gene.func.gii')).agg_data().astype(float)
    vmin, vmax = np.nanpercentile(lh, 2), np.nanpercentile(lh, 98)
    for view in ['lateral', 'medial']:
        fig = plt.figure(figsize=(2.4, 2.0))
        ax = fig.add_subplot(111, projection='3d')
        plotting.plot_surf(infl_l, lh, hemi='left', view=view, cmap=cmap,
                           vmin=vmin, vmax=vmax, colorbar=False,
                           axes=ax, figure=fig)
        fig.savefig(os.path.join(tmp_dir, f'{g}_{view}.png'), dpi=150,
                    bbox_inches='tight', facecolor='white')
        plt.close(fig)
    print('rendered', g, flush=True)

# ---- stage 2: composite grid with imshow (fast) ----
import matplotlib.image as mpimg
n = len(genes)
fig, axes = plt.subplots(n, 2, figsize=(5.0, 1.5 * n))
fig.patch.set_facecolor('white')
for i, g in enumerate(genes):
    for j, view in enumerate(['lateral', 'medial']):
        ax = axes[i, j]
        ax.imshow(mpimg.imread(os.path.join(tmp_dir, f'{g}_{view}.png')))
        ax.axis('off')
    axes[i, 0].text(-0.08, 0.5, gene_labels[g], transform=axes[i, 0].transAxes,
                    ha='right', va='center', fontsize=13, fontweight='bold')

cbar_ax = fig.add_axes([0.93, 0.25, 0.015, 0.5])
sm = cm.ScalarMappable(cmap=cmap); sm.set_array([0, 1])
cb = fig.colorbar(sm, cax=cbar_ax, ticks=[0, 1])
cb.ax.set_yticklabels(['low', 'high']); cb.set_label('AHBA expression', fontsize=10)

fig.suptitle('hormone / ED receptor gene expression', fontsize=15, fontweight='bold', y=0.995)
fig.subplots_adjust(left=0.1, right=0.9, top=0.96, bottom=0.01, wspace=-0.3, hspace=0.05)
out = os.path.join(out_dir, 'HORMONEgenes_surface_expression.png')
fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print('Saved', out, flush=True)
