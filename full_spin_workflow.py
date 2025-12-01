import neuromaps
from neuromaps import nulls, datasets, images
from neuromaps.stats import compare_images
from neuromaps.images import load_data
from neuromaps import resampling
from nilearn import plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import visualizations_functions as vf

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
source_map = {'source':'hill2010', 'desc':'evoexp', 'space':'fsLR', 'den':'164k'}

# Full implementation of the Alexander Bloch Spin Test
def full_spin_test(src: dict, trg: dict):
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
    rotated = nulls.alexander_bloch(
        src_sphere,
        atlas = trg_space,
        density = trg_den,
        n_perm = 1000,
        seed = 1234
    )
    r, p, null = compare_images(
        src_sphere,
        trg_sphere,
        metric='pearsonr',
        nulls=rotated,
        return_nulls=True
    )
    end = time.perf_counter()
    time_elapsed = end - start
    formal_name = vf.map_names.get(trg_title)
    results_dict = {'target map':formal_name, 
                'r_emp':r, 
                'p_spin':p, 
                'nulls':null,
                'runtime': time_elapsed
               }
    return results_dict

if __name__ == "__main__":
    results = []
    for target_map in target_maps:
        results.append(full_spin_test(source_map, target_map))
    vf.make_multiTest_table(pd.DataFrame(results), alpha=0.05, p_col='p_spin', label_col='target map')
    vf.plot_box_plot(results)

    # # uncomment this section to plot brain maps
    # for trg_map in target_maps:
    #     vf.plot_brain_map(trg_map, vf.map_names, vf.brain_map_settings)
