import neuromaps
from neuromaps import datasets
from neuromaps.images import load_data
from neuromaps.datasets import fetch_atlas
from nilearn import plotting
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import time
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

map_names = {
    'genepc1': 'PC1 Gene Expression',
    'myelinmap': 'T1w/T2w Ratio',
    'thickness': 'Cortical Thickness',
    'devexp': 'Developmental Expansion',
    'fcgradient01': 'Functional Gradient',
    'intersubjvar': 'Intersubject Variability',
    'cbf': 'Cerebral Blood Flow',
    'cbv': 'Cerebral Blood Volume',
    'cmr02': 'Oxygen Metabolism',
    'cmrglc': 'Glucose Metabolism',
    'scalingnih': 'Allometric Scaling (NIH)',
    'scalingpnc': 'Allometric Scaling (PNC)',
    'evoexp': 'Evolutionary Expansion'
}
brain_map_settings = {
    'evoexp': {'cmap': 'inferno', 'vmin': -2.7, 'vmax': 2.7},

    'genepc1': {'cmap': 'magma', 'vmin': -2.7, 'vmax': 2.7},
    'myelinmap':{'cmap': 'viridis', 'vmin': None, 'vmax': None},
    'thickness': {'cmap': 'viridis', 'vmin': None, 'vmax': None},
    'devexp': {'cmap': 'inferno', 'vmin': -1, 'vmax': 1},
    'fcgradient01': {'cmap': 'rainbow', 'vmin': None, 'vmax': None},
    'intersubjvar': {'cmap': 'inferno', 'vmin': None, 'vmax': None},
    
    'cbf': {'cmap': 'viridis', 'vmin': 'special_perc', 'vmax': 'special_perc'},
    'cbv': {'cmap': 'viridis', 'vmin': None, 'vmax': None},
    'cmr02': {'cmap': 'viridis', 'vmin': 'special_perc', 'vmax': 'special_perc'},
    'cmrglc': {'cmap': 'viridis', 'vmin': 'special_perc', 'vmax': 'special_perc'},
    'scalingnih': {'cmap': 'seismic', 'vmin': None, 'vmax': None},
    'scalingpnc': {'cmap': 'seismic', 'vmin': None, 'vmax': None},
}

## Box plot function
def plot_box_plot(results):
    results_df = pd.DataFrame(results)
    box_data = [np.ravel(np.array(n)) for n in results_df['nulls']]
    positions = np.arange(1, len(box_data) + 1)

    fig, ax = plt.subplots(figsize=(12, 8))
    #boxplots for target maps null distributions
    bp = ax.boxplot(
        box_data,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        boxprops=dict(facecolor='white', edgecolor='black', linewidth=1),
        medianprops=dict(color='black', linewidth=1),
        whiskerprops=dict(color='black'),
        capprops=dict(color='black'),
        flierprops=dict(marker='o', color='gray', markersize=3, alpha=0.4)
    )

    #add spin test r correlation values
    for i, (r, p) in enumerate(zip(results_df['r_emp'], results_df['p_spin'])):
        color = 'red' if p < 0.05 else '#e6a67a'  # red = significant, orange = non-significant
        ax.scatter(
            positions[i],
            r,
            color=color,
            s=80,
            edgecolor='black',
            linewidth=0.5,
            zorder=3
        )

    ax.axhline(0, color='gray', linestyle='--', linewidth=1)
    ax.set_xticks(positions)
    ax.set_xticklabels(results_df['target map'], rotation=40, ha='right', fontsize=10)
    ax.set_ylabel("Pearson's r", fontsize=12)
    ax.set_xlabel("Target maps", fontsize=12)
    ax.set_ylim(-0.8, 0.8)
    ax.set_title("Empirical correlations vs spatial nulls", fontsize=13, pad=15)

    #legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w',
           label=r'Empirical ($P_{spin} ≥ 0.05$)',
           markerfacecolor='#e6a67a', markeredgecolor='black', markersize=8),
        Line2D([0], [0], marker='o', color='w',
           label=r'Empirical ($P_{spin} < 0.05$)',
           markerfacecolor='red', markeredgecolor='black', markersize=8),
        Line2D([0], [0], color='black', lw=1, label='Spatial null')
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=False)
    plt.show()

## Function to make summary table for multiple comparisons
def make_multiTest_table(results_df, alpha=0.05, p_col='p_spin', label_col='target map'):

    p_vals = np.array(results_df[p_col])

    reject_fdr_bh, p_fdr_bh, _, _ = multipletests(p_vals, alpha=alpha, method='fdr_bh')
    reject_fdr_by, p_fdr_by, _, _ = multipletests(p_vals, alpha=alpha, method='fdr_by')
    reject_bonf, p_bonf, _, _ = multipletests(p_vals, alpha=alpha, method='bonferroni')

    # Add to original df (if you want to keep them there too)
    results_df = results_df.copy()
    results_df['p_Bonf'] = p_bonf
    results_df['reject_Bonf'] = reject_bonf
    results_df['p_FDR_BH'] = p_fdr_bh
    results_df['reject_FDR_BH'] = reject_fdr_bh
    results_df['p_FDR_BY'] = p_fdr_by
    results_df['reject_FDR_BY'] = reject_fdr_by
    results_df['reject_uncorrected'] = results_df[p_col] < alpha

    # Build summary table in your requested column order
    summary_df = results_df[[ 
        label_col,
        p_col,
        'reject_uncorrected',
        'p_Bonf',
        'reject_Bonf',
        'p_FDR_BH',
        'reject_FDR_BH',
        'p_FDR_BY',
        'reject_FDR_BY'
    ]].copy()

    summary_df.columns = [
        'Target map',
        'p (raw)',
        'Significant (p<0.05)',
        'p_Bonferroni (corrected)',
        'Significant (Bonf<0.05)',
        'p_FDR (BH corrected)',
        'Significant (FDR BH<0.05)',
        'p_FDR (BY corrected)',
        'Significant (FDR BY<0.05)'
    ]

    summary_df = summary_df.round(4)
    return summary_df

# Function to creat brain maps
def plot_brain_map(map: dict, map_names: dict, brain_map_settings: dict):
    """
    Plots the brain map given in the map dictionary
    Parameters:
    map: the map as a dictionary with the needed parameters for feth_annotation
    map_names: dictionary of formal names for readability and plot titles
    Outputs: A brain map plot
    """
    start_time  = time.perf_counter() #timing just for user info

    map_paper, map_desc, map_space, map_den = map.values()
    #fetch source map and target map files
    src_map = datasets.fetch_annotation(**map)

    settings = brain_map_settings.get(map_desc, {})
    cmap = settings.get('cmap', 'inferno')

    fig = plt.figure(figsize=(10, 4))
    fslr = fetch_atlas(map_space, map_den)
    surf_mesh_left = fslr['inflated'].L
    surf_mesh_right = fslr['inflated'].R
    data_full = load_data(src_map)

    if settings.get('vmin') ==  'special_perc' and settings.get('vmax') == 'special_perc':
        vmin, vmax = np.percentile(data_full[~np.isnan(data_full)], [10, 95])
    elif settings.get('vmin') is not None and settings.get('vmax') is not None:
        vmin, vmax = settings['vmin'], settings['vmax']
    else:
        vmin, vmax = np.percentile(data_full[~np.isnan(data_full)], [3, 98])

    #if target map has both hemispheres, plot both hemispheres lateral
    if(len(src_map)==2):
        data_l = load_data(src_map[0])
        data_r = load_data(src_map[1])
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        plotting.plot_surf(
            surf_mesh=surf_mesh_left,
            surf_map=data_l,
            hemi='left',
            view='lateral',
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            colorbar=False,
            axes=ax1,
            title='Left hemisphere'
        )
        # right hemi
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        plotting.plot_surf(
            surf_mesh=surf_mesh_right,
            surf_map=data_r,
            hemi='right',
            view='lateral',
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            colorbar=False,
            axes=ax2,
            title='Right hemisphere'
        )
        # color bar
        sm = plt.cm.ScalarMappable(cmap=cmap)
        sm.set_clim(vmin, vmax)
        cbar = fig.colorbar(sm, ax=[ax1, ax2], shrink=0.6, location='right')
        cbar.set_label(f"{map_names.get(map_desc)}({map_space} {map_den})", fontsize=11)
        plt.suptitle(f"{map_names.get(map_desc)}", fontsize=14)
        plt.show()


    else: #if only one hemisphere, only plot that hemisphere lateral and medial
        plot_kwargs = dict(
            surf_mesh=surf_mesh_right,
            surf_map=data_full,
            hemi='right',
            bg_on_data=True,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            colorbar=False
        )
        #plot right veiw
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        plotting.plot_surf(
            view='lateral',
            axes=ax1,
            **plot_kwargs
        )
        #plot from left view
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        plotting.plot_surf(
            view='medial',
            axes=ax2,
            **plot_kwargs
        )
        #color bar
        sm = plt.cm.ScalarMappable(cmap=cmap)
        sm.set_clim(vmin, vmax)
        cbar = fig.colorbar(sm, ax=[ax1, ax2], shrink=0.6, location='right')
        cbar.set_label(f"{map_names.get(map_desc)}({map_space} {map_den})", fontsize=11)
        plt.suptitle(f"{map_names.get(map_desc)}", fontsize=14)
        plt.show()

    end_time = time.perf_counter()
    print(f"Total time for plotting: {(end_time - start_time):.2f} seconds")

    return None

