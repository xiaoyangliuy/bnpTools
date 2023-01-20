#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 06:29:29 2020

@author: yluo89
"""


import os, openslide, collections
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

def loadNDPIBestResSingleFile(fname):
    lv, dsmp, region = None, None, None
    with openslide.OpenSlide(fname) as slide:
        level = slide.level_count
        dim = slide.level_dimensions
        ds = slide.level_downsamples
        idx = np.argmin(ds)
        region = slide.read_region((0,0),idx,dim[idx])
        lv = idx
        dsmp = ds[idx]
    return lv, dsmp, region

def loadNDPIDir(fdir):
    data = collections.defaultdict(list)
    fnames = os.listdir(fdir)
    for f in fnames:
        if '.ndpi' in f:
            lv, dsmp, region = loadNDPIBestResSingleFile(os.path.join(fdir,f))
            data['name'].append(f[:-5])
            data['level'].append(lv)
            data['downsample'].append(dsmp)
            data['img'].append(region)
    return data


