from neuromaps import nulls, datasets, images
from neuromaps.stats import compare_images
from neuromaps import resampling
from nilearn import plotting
import matplotlib.pyplot as plt
import numpy as np
import time
target_maps = [
    {'source':'raichle', 'desc':'cbf', 'space':'fsLR', 'den':'164k'},
    {'source':'raichle', 'desc':'cbv', 'space':'fsLR', 'den':'164k'},
    {'source':'raichle', 'desc':'cmr02', 'space':'fsLR', 'den':'164k'},
    {'source':'raichle', 'desc':'cmrglc', 'space':'fsLR', 'den':'164k'},
    {'source':'reardon2018', 'desc':'scalingnih', 'space':'civet', 'den':'41k'},
    {'source':'reardon2018', 'desc':'scalingpnc', 'space':'civet', 'den':'41k'}
]
source_map = {'source':'hill2010', 'desc':'evoexp', 'space':'fsLR', 'den':'164k'}
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

results = []
for target_map in target_maps:
    results.append(full_spin_test(source_map, target_map))

