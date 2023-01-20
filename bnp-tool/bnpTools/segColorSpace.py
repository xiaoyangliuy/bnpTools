#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 09:45:57 2020

@author: yluo89
Separate color space of HED (hematoxylin, eosin, and DAB) on stained tissue
"""

import os, cv2
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors


# hc and lc are the upper and lower limit for segmentation
# in the case of using hc or lc for hue value, the expected inputs are (h_hue, 255,255) and (l_hue, 0, 0)
def colorSegamentation(roi, hc1, lc1, hc2, lc2, showPlots=True, figsize=(3,5)):
    roihsv = cv2.cvtColor(roi, cv2.COLOR_RGB2HSV)
    c1mask = getMask(roihsv, lc1, hc1)
    c2mask = getMask(roihsv, lc2, hc2)
    cmap1 = (255-cv2.cvtColor(roi, cv2.COLOR_RGB2GRAY)) * c1mask
    cmap2 = (255-cv2.cvtColor(roi, cv2.COLOR_RGB2GRAY)) * c2mask
    
    if showPlots is True:
        pltstr=['roi','cmap1','cmap2',['cmap1','cmap2']]
        cmap=[None,'Blues', 'Reds', ['Blues','Reds']]
        fig,axes = plt.subplots(2,2,figsize=figsize)
        for s_, ax_, c_ in zip(pltstr, axes.ravel(), cmap):
            if len(s_) != 2:
                d = eval(s_)
                ax_.imshow(d, cmap=c_)
            else:
                d1 = eval(s_[0])
                d2 = eval(s_[1])
                ax_.imshow(d1, cmap=c_[0])
                ax_.imshow(d2, cmap=c_[1], alpha=0.4)
            ax_.axis('off')
        plt.tight_layout()
    return c1mask, c2mask, cmap1, cmap2, fig


def hsvDistribution(roi, elv = None, angle=None):
    if angle is None:
        angle = [45,0,90]
        
    if elv is None:
        elv = 30
        
    roihsv = cv2.cvtColor(roi, cv2.COLOR_RGB2HSV)
    fig = plt.figure(figsize=(10,10))
    for i, a in enumerate(angle):
        ax = fig.add_subplot(1, len(angle),i+1, projection='3d')
        splitColorSpace(roi, roihsv, labels=['Hue','Saturation','Value'], fig=fig, 
                    axis=ax, elv=elv, angle=a)
    
    plt.tight_layout()
    return fig

def rgbDistribution(roi, elv = None, angle=None):
    if angle is None:
        angle = [45,0,90]
        
    if elv is None:
        elv = 30
        
    fig = plt.figure(figsize=(10,10))
    for i, a in enumerate(angle):
        ax = fig.add_subplot(1, len(angle),i+1, projection='3d')
        splitColorSpace(roi, roi, labels=['Red','Green','Blue'], fig=fig, 
                    axis=ax, elv=elv, angle=a)
    
    plt.tight_layout()
    return fig

def splitColorSpace(img,img_hsv,labels = None, axis=None, fig=None, elv=None, angle=None):
    if labels is None:
        labels = ['a', 'b', 'c']
        
    a, b, c = cv2.split(img_hsv)
    if (axis is None) | (fig is None):
        fig = plt.figure()
        axis = fig.add_subplot(1, 1, 1, projection="3d")
        
    pixel_colors = img.reshape((img.shape[0]*img.shape[1], 3))
    norm = colors.Normalize(vmin=-1.,vmax=1.)
    norm.autoscale(pixel_colors)
    pixel_colors = norm(pixel_colors).tolist()
    axis.scatter(a.flatten(), b.flatten(), c.flatten(), facecolors=pixel_colors, marker=".")
    axis.set_xlabel(labels[0])
    axis.set_ylabel(labels[1])
    axis.set_zlabel(labels[2])
    rotate3Dplot(fig, axis, elv=elv, angle=angle)
    return fig,axis
    
def rotate3Dplot(fig,axis,elv=None, angle=None):
    if elv is None: elv = 0
    if angle is None: angle = 0
    axis.view_init(elv,angle)
    fig
    return fig
    
def getMask(roihsv,lim1, lim2):
    mask = cv2.inRange(roihsv, lim1, lim2)
    mask[mask==255] = 1
    return mask

def plotmaps(roi, result, figsize=None, cmap=None, axisOff=True):
    # r,g,b = cv2.split(roi)
    fig,axes = plt.subplots(1,2,figsize=figsize)
    axes[0].imshow(roi)
    axes[1].imshow(result,cmap=cmap)
    if axisOff is True:
        axes[0].axis('off')
        axes[1].axis('off')
    return fig