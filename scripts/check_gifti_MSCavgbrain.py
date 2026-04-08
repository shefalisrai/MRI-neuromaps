"""
Plots.func.gii z-stat maps on the MSC inflated surface template that I have.
SR
"""

import os
os.environ['PATH'] += ':/Applications/workbench/bin_macosx64'
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nilearn import plotting
from neuromaps import transforms

lh_surf = '/Users/shefalirai/Desktop/Docs/sub-MSCAvg.L.inflated.32k_fs_LR.surf.gii'
rh_surf = '/Users/shefalirai/Desktop/Docs/sub-MSCAvg.R.inflated.32k_fs_LR.surf.gii'

gifti_dir = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/3Groupcontrasts_funcgiftis/'
out_dir      = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/surface_checks/'
resamp_dir   = os.path.join(out_dir, 'resampled_32k')
os.makedirs(out_dir,    exist_ok=True)
os.makedirs(resamp_dir, exist_ok=True)

measures =['T1w_contrast', 'sulc', 'thickness', 'area','RND-gm','RNI-gm', 'RND-wm','RNI-wm']
contrasts = ['EDR_HC', 'EDBP_HC', 'EDR_EDBP']

def inspect_gifti(fpath):
    img  = nib.load(fpath)
    data = img.darrays[0].data.astype(np.float32)
    print(f'    {os.path.basename(fpath)}')
    print(f'      vertices : {len(data)}')
    print(f'      range    : [{np.nanmin(data):.3f}, {np.nanmax(data):.3f}]')
    print(f'      NaNs     : {np.sum(np.isnan(data))}')
    print(f'      zeros    : {np.sum(data == 0)} ({100*np.mean(data == 0):.1f} %)')
    return data


def resample_to_fslr32k(lh_gii, rh_gii, tag):
    """
    Resample a (lh, rh) pair from fsaverage ico5 to fs_LR 32k using neuromaps.
    Returns paths to the resampled files 
    """
    lh_out = os.path.join(resamp_dir, f'lh_{tag}_32k.func.gii')
    rh_out = os.path.join(resamp_dir, f'rh_{tag}_32k.func.gii')

    if os.path.exists(lh_out) and os.path.exists(rh_out):
        print('  [resample] cached, skipping')
        return lh_out, rh_out

    print('  [resample] fsaverage ico5 -> fs_LR 32k ...')

    lh_r, rh_r = transforms.fsaverage_to_fslr(
        (lh_gii, rh_gii),
        target_density='32k',
    )

    # lh_r / rh_r are nibabel GiftiImage objects 
    nib.save(lh_r, lh_out)
    nib.save(rh_r, rh_out)
    print(f'    saved resampled LH -> {os.path.basename(lh_out)}')
    print(f'    saved resampled RH -> {os.path.basename(rh_out)}')
    return lh_out, rh_out


def plot_four_views(lh_gii, rh_gii, title, out_png):
    """2x2 figure: LH lateral, LH medial, RH lateral, RH medial."""
    lh_data = nib.load(lh_gii).darrays[0].data.astype(np.float32)
    rh_data = nib.load(rh_gii).darrays[0].data.astype(np.float32)

    vmax = float(np.percentile(np.abs(np.concatenate([lh_data, rh_data])), 98))
    vmax = max(vmax, 2.0)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9),
                             subplot_kw={'projection': '3d'})
    fig.patch.set_facecolor('white')

    panels = [
        ('lateral', 'left',  lh_gii, lh_surf, axes[0, 0]),
        ('medial',  'left',  lh_gii, lh_surf, axes[0, 1]),
        ('lateral', 'right', rh_gii, rh_surf, axes[1, 0]),
        ('medial',  'right', rh_gii, rh_surf, axes[1, 1]),
    ]

    for view, hemi, gii_f, surf_f, ax in panels:
        plotting.plot_surf_stat_map(
            surf_mesh=surf_f,
            stat_map=gii_f,
            hemi=hemi,
            view=view,
            bg_map=None,           # no sulc 
            colorbar=True,
            vmax=vmax,
            cmap='RdBu_r',
            symmetric_cbar=True,
            axes=ax,
            title=f'{hemi} {view}',
        )

    fig.suptitle(title, fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  -> saved: {out_png}')

print('Checking template surfaces...')
for surf_path in [lh_surf, rh_surf]:
    if not os.path.exists(surf_path):
        raise FileNotFoundError(f'Template not found: {surf_path}')
    coords, _ = nib.load(surf_path).agg_data()
    print(f'  {os.path.basename(surf_path)}  --  {len(coords)} vertices  OK')
print()

for measure in measures:
    for contrast in contrasts:

        lh_gii = os.path.join(
            gifti_dir,
            f'lh.postpsubs_vertexwise_{contrast}_{measure}_sm256_baseline_3Group.zstat.func.gii'
        )
        rh_gii = os.path.join(
            gifti_dir,
            f'rh.postpsubs_vertexwise_{contrast}_{measure}_sm256_baseline_3Group.zstat.func.gii'
        )

        if not os.path.exists(lh_gii) or not os.path.exists(rh_gii):
            print(f'[skip] missing files for {measure} -- {contrast}')
            continue

        print(f'\n{measure} | {contrast}')
        print('  Original (fsaverage ico5):')
        inspect_gifti(lh_gii)
        inspect_gifti(rh_gii)

        tag = f'{measure}_{contrast}'
        lh_32k, rh_32k = resample_to_fslr32k(lh_gii, rh_gii, tag)

        print('  Resampled (fs_LR 32k):')
        inspect_gifti(lh_32k)
        inspect_gifti(rh_32k)

        out_png = os.path.join(out_dir, f'{measure}_{contrast}.png')
        plot_four_views(lh_32k, rh_32k,
                        title=f'{measure}  |  {contrast}',
                        out_png=out_png)

print('\nAll done. PNGs in:', out_dir)
