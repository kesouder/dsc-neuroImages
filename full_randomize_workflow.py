from neuromaps import nulls, datasets, images
from neuromaps.stats import compare_images
from neuromaps import resampling
from nilearn import plotting
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import time
import pandas as pd
target_maps = [
    {'source':'abagen', 'desc':'genepc1', 'space':'fsaverage', 'den':'10k'},
    {'source':'hcps1200', 'desc':'myelinmap', 'space':'fsLR', 'den':'32k'},
    {'source':'hcps1200', 'desc':'thickness', 'space':'fsLR', 'den':'32k'},
    {'source':'hill2010', 'desc':'devexp', 'space':'fsLR', 'den':'164k'},
    {'source':'margulies2016', 'desc':'fcgradient01', 'space':'fsLR', 'den':'32k'},
    {'source':'mueller2013', 'desc':'intersubjvar', 'space':'fsLR', 'den':'164k'},
    {'source':'raichle', 'desc':'cbf', 'space':'fsLR', 'den':'164k'},
    {'source':'raichle', 'desc':'cbv', 'space':'fsLR', 'den':'164k'},
    {'source':'raichle', 'desc':'cmr02', 'space':'fsLR', 'den':'164k'},
    {'source':'raichle', 'desc':'cmrglc', 'space':'fsLR', 'den':'164k'},
    {'source':'reardon2018', 'desc':'scalingnih', 'space':'civet', 'den':'41k'},
    {'source':'reardon2018', 'desc':'scalingpnc', 'space':'civet', 'den':'41k'}
]
# Official names for plotting later
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
source_map = {'source':'hill2010', 'desc':'evoexp', 'space':'fsLR', 'den':'164k'}
def full_random_test(src: dict, trg: dict):
    src_paper, src_title, src_space, src_den = src.values()
    trg_paper, trg_title, trg_space, trg_den = trg.values()
    #fetch source map and target map files
    start = time.perf_counter()
    src_map = datasets.fetch_annotation(**src)
    trg_map = datasets.fetch_annotation(**trg)
    #if target map have both hemispheres, use the right one
    if(len(trg_map)==2):
        trg_map = trg_map[1]
    if src_den != trg_den:
        src_res, trg_res = resampling.resample_images(
            src_map,
            trg_map,
            src_space=src_space,
            trg_space=trg_space,
            hemi='R',
            resampling='downsample_only' #sample to size of target density
        )
        src_data = images.load_data(src_res)
        trg_data = images.load_data(trg_res)
    else:
        src_data = images.load_data(src_map)
        trg_data = images.load_data(trg_map)
    #Create nan values for left brain
    L_nan = np.full_like(src_data, np.nan)
    src_sphere = np.hstack([L_nan, src_data])
    trg_sphere = np.hstack([L_nan, trg_data])
    randomized = []
    for i in range(1000):
        random_shuffle_R = src_data[np.random.permutation(len(src_data))]
        randomized.append(np.hstack([L_nan, random_shuffle_R]))
    randomized = np.array(randomized).T
    r, p, null = compare_images(
        src_sphere,
        trg_sphere,
        metric='pearsonr',
        nulls=randomized,
        nan_policy='omit',
        return_nulls=True
    )
    end = time.perf_counter()
    time_elapsed = end - start
    results_dict = {'target map':trg_title, 
                'r_emp':r, 
                'p_spin':p, 
                'nulls':null,
                'runtime': time_elapsed
               }
    return results_dict

def plot_box_plot(results):
    results_df = pd.DataFrame(results)
    box_data = [np.ravel(np.array(n)) for n in results_df['nulls']]
    positions = np.arange(1, len(box_data) + 1)

    fig, ax = plt.subplots(figsize=(11, 8))
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
    plt.savefig('Randomized Permutation Box Plot.png')
    plt.show()


if __name__ == "__main__":
    results = []
    for target_map in target_maps:
        results.append(full_random_test(source_map, target_map))
    plot_box_plot(results)

