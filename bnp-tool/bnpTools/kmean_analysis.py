#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:55:21 2020

@author: yluo89
"""


import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
import skimage.filters
import os, cv2, h5py
from matplotlib.patches import Rectangle



# K-mean cluster analysis
def kmean_analysis(n_clusters, data, random_state, sigma = None, cval = None,
                   plotoption = None, savefig = None, fname = None):
    
    data[np.isnan(data)] = 1e-5
    data[np.isinf(data)] = 1e-5
    
    if sigma is None: sigma = 1
    data_blur = skimage.filters.gaussian(data, sigma = sigma)
    
    km = KMeans(n_clusters = n_clusters,random_state=random_state)
    km.fit(data_blur.reshape(-1,1))
    
    km_label = np.reshape(km.labels_, data.shape)
    
    # sort label based on center
    srtIndex = np.argsort(km.cluster_centers_[:,0])
    for i, s in enumerate(srtIndex):
        km_label[km_label == s] = -(i+1)
    km_label = np.multiply(-1,km_label)-1
    km_bool = km_label.copy()
    km_bool[km_bool > 1] = 1
    
    fig = None
    
    if plotoption:
        fig, ax = plt.subplots(1,4,figsize=(7,4))
        a = ax[0].imshow(data, cmap = plt.cm.get_cmap('Greys_r'))
        if cval == None:
            cval = a.get_clim()
        else:
            a.set_clim(cval)
        c = ax[1].imshow(data_blur, cmap = plt.cm.get_cmap('inferno'))
        k = ax[2].imshow(km_label, vmin = 0, vmax = n_clusters-1)
        b = ax[3].imshow(np.multiply(data,km_bool), cmap = plt.cm.get_cmap('Greys_r'))
        b.set_clim(cval)
        c.set_clim(cval)
        fig.colorbar(a, ax=ax[0], orientation='horizontal', shrink = 0.8)
        fig.colorbar(c, ax = ax[1], orientation='horizontal', shrink = 0.8)
        fig.colorbar(k, ax = ax[2], orientation='horizontal', shrink = 0.8)
        fig.colorbar(b, ax = ax[3], orientation = 'horizontal', shrink = 0.8)
        
        
        map_label = ['data','blur', 'blur-kmean', 'data * kmean']
        for ax_, l in zip(ax, map_label):
            ax_.axis('off')
            ax_.axis('equal')
            ax_.set_title(l)
            # ax_.axis('scaled')
#            ax_.text(2, 7, l, color = 'w')
        plt.tight_layout()
        plt.show()
        
        if (savefig == 1) & (fname is not None):
            fig.savefig(fname, dpi = 300)
    
    return km, km_label, km_bool, fig


