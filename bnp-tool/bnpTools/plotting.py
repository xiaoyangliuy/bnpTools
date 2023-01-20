#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 08:11:52 2020

@author: yluo89
"""
import matplotlib.pyplot as plt
import numpy as np
# from sklearn.cluster import KMeans
import skimage.filters
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
# import cv2

def plotProj(bnp_data, elm, figsize=(25, 22), use_imshow=False, ncol = 8, vmax = None, n_std = 2, colorbar = False, cmap = 'RdYlBu_r'):
    nrow = int(np.ceil(len(bnp_data)/ncol))
    fig,ax = plt.subplots(nrow, ncol, figsize=figsize) 
    axes = ax.ravel()
    counter = 0
    clim = []
    for i, row in bnp_data.iterrows():
        dat = row[elm]
        dat[np.isinf(dat) | np.isnan(dat)] = 0
        
        if vmax is None:
            vmax1 = np.nanmean(dat) + n_std * np.nanstd(dat)
        else:
            vmax1 = vmax
        
        if use_imshow:
            img = axes[counter].imshow(dat,vmin = 0,
                       vmax = vmax1,
                       cmap = cmap)
        else:
            img = axes[counter].pcolormesh(row['xval'], row['yval'], dat,vmin = 0,
                                   vmax = vmax1,
                                   cmap = cmap, shading='auto')

        axes[counter].set_title('%.2f (ss: %.2f) \n %s'%(row['theta'], np.abs(np.diff(row['xval']))[0], row['scan_num']), 
                                fontsize=18)
        axes[counter].set_xticks([])
        axes[counter].set_yticks([])
        
        if colorbar:
            divider = make_axes_locatable(axes[counter])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = fig.colorbar(img, cax = cax, shrink=0.8)
            cbar.ax.tick_params(labelsize=14)
            
        clim.append(img.get_clim())
        counter += 1
        
    for i in range(len(axes)-counter):
        axes[i+counter].axis('off')
    print('Total number of useful projections: %d'%(counter))
    plt.tight_layout()
    return fig, ax, clim

def addColorBar(fig, img, ax):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(img, cax = cax, shrink=0.8)
    cbar.ax.tick_params(labelsize=12)
    return cbar

def plotProjArray(bnp_data, elm_idx, angles, stepsize = None, figsize=(25, 22), ncol = 8, vmax = None, colorbar = False, cmap = 'inferno'):
    nrow = int(np.ceil(bnp_data.shape[2]/ncol))
    fig,ax = plt.subplots(nrow, ncol, figsize=figsize) 
    axes = ax.ravel()
    counter = 0
    clim = []
    for i in range(bnp_data.shape[2]):
        dat = bnp_data[elm_idx,:,i,:]
        dat[np.isinf(dat) | np.isnan(dat)] = 0
        if vmax is None:
            vmax1 = np.nanmean(dat) + 2 * np.nanstd(dat)
        else:
            vmax1 = vmax
        
        img = axes[counter].imshow(dat,vmin = 0,
                   vmax = vmax1,
                   cmap = cmap, origin='lower')

        if stepsize is None: ss = np.nan
        else: ss = float(stepsize[i])
            
        if colorbar:
            divider = make_axes_locatable(axes[counter])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = fig.colorbar(img, cax = cax, shrink=0.8)
            cbar.ax.tick_params(labelsize=14)
            
        axes[counter].set_title('%.2f (ss: %.2f)'%(angles[i], ss), fontsize=18)
        axes[counter].set_xticks([])
        axes[counter].set_yticks([])
        clim.append(img.get_clim())
        counter += 1
        
    for i in range(len(axes)-counter):
        axes[i+counter].axis('off')
    plt.tight_layout()
    return fig, ax, clim

# Overlay two maps
def xrfOverlay(map1, map2, vmax_m1 = None, vmax_m2 = None, labels = None, ax = None, cmap=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.gca()
    else:
        fig = 0
    
    if cmap is None:
        cmap = ['Reds','Greens']
    
    ax.imshow(map1, cmap=cmap[0], alpha = 1, vmax = vmax_m1)
    ax.imshow(map2, cmap =cmap[1], alpha = 0.5, vmax = vmax_m2)
    
    if labels is not None:
        ax.text(map2.shape[1]*0.05, map2.shape[0]*0.85, labels[0], color=plt.get_cmap(cmap[0])(0.8), fontsize=12)
        ax.text(map2.shape[1]*0.05, map2.shape[0]*0.95, labels[1], color=plt.get_cmap(cmap[1])(0.8), fontsize=12)
        
    ax.set_xticks([])
    ax.set_yticks([])
    
    return fig

def plotMeanIntensityTheta(bnp_data, img_tags, img_def, elm_chs, ax, 
                           cmap=None, cval = None, savefig = None, 
                           fname = None):
    
    if cmap == None: cmap = plt.cm.get_cmap('RdYlBu')
    if cval == None: cval = np.arange(0.1, 0.8, 0.7/(len(elm_chs)))
    
    for i, img_tag in enumerate(img_tags):
        data_pd = bnp_data[img_tag]
        for e, ax_ in zip(elm_chs, ax.ravel()):
            xval = data_pd['theta']
            yval = np.array([np.nanmean(i) for i in data_pd[e]])
            ax_.plot(xval, yval, label=img_def[i], marker='.', 
                     linewidth=0, markersize=5,
                     color=cmap(cval[i]),alpha = 0.8)
            ax_.set_yscale('log')
            ax_.set_title(e)

    ax_.legend(loc='lower left', bbox_to_anchor=(1, 0.15), 
               frameon=False, markerscale=2)
    
    if len(ax.shape) == 1:
        ax[0].set_ylabel('Mean Intensity')
        for i in range(ax.shape[0]):
            ax[i].set_xlabel('Theta ($^\circ$)')
    else:
        for i in range(ax.shape[0]):
            ax[i,0].set_ylabel('Mean Intensity')
        for i in range(ax.shape[1]):
            ax[-1,i].set_xlabel('Theta ($^\circ$)')
        
    plt.tight_layout()
    plt.show()
    
    if (savefig == 1) & (fname is not None):
        plt.savefig(fname, dpi = 300)
    
    
def plotMapDiff(map1, map2, theta1=None, theta2=None, vmax = None, 
                savefig = None, fname = None, ax = None):    
    
    cmap = plt.cm.get_cmap('coolwarm')
    xloc = map1.shape[1]/2+map1.shape[1]*0.15
    yloc = map1.shape[0]-map1.shape[0]*0.1
    
    if ax is None:
        fig = plt.figure()
        plt.imshow(map1 * (-1), cmap = cmap, 
                   alpha = 0.5, vmax = vmax)
        plt.imshow(map2, cmap = cmap, 
                   alpha = 0.5, vmax = vmax)
        plt.axis('off')
        if (theta1 is not None) & (theta2 is not None):
            plt.text(xloc, yloc, '%.2f'%(theta1), color = cmap(0.3))
            plt.text(xloc, yloc+map1.shape[0]*0.05, '%.2f'%(theta2), color = cmap(0.8))
        if savefig:
            plt.savefig(fname, dpi = 300)
        
    else:
        ax.imshow(map1 * (-1), cmap = cmap, 
                   alpha = 0.5, vmax = vmax)
        ax.imshow(map2, cmap = cmap, 
                   alpha = 0.5, vmax = vmax)
        ax.axis('off')
        if (theta1 is not None) & (theta2 is not None):
            ax.text(xloc, yloc, '%.2f'%(theta1), color = cmap(0.3))
            ax.text(xloc, yloc+map1.shape[0]*0.05, '%.2f'%(theta2), color = cmap(0.8))
    
    
def plot2map(data1, data2, savefig=None, fname = None):
    fig, axes = plt.subplots(1,2,figsize=(5,2.5))
    a = axes[0].imshow(data1, cmap = plt.cm.get_cmap('Greys_r'))
    b = axes[1].imshow(data2, cmap = plt.cm.get_cmap('Greys_r'))
    cval = a.get_clim()
    b.set_clim(cval)
    
    if (savefig) & (fname is not None):
        plt.tight_layout()
        fig.savefig(fname, dpi = 300)
    
    return fig, axes

def plotSingleChMaps(df, ch_name):
    for i, row in df.iterrows():
        plt.figure()
        plt.imshow(row[ch_name])
        
def plot4DimProj(proj_data, index):
    counter = 0
    for i in proj_data[index,:,:,:]:
        plt.figure()
        plt.imshow(i)
        plt.title(str(counter))
        counter += 1
        
def plot3DimProj(proj_data):
    counter = 0
    for i in range(len(proj_data)):
        d = proj_data[i,:,:]
        plt.figure()
        plt.imshow(d)
        plt.title(str(counter))
        counter += 1
        
def plotRecon(recon, index=None, savefig = None, outdir = None, smp = None, 
              interpolation = None):
    a,b,c = recon.shape
    numplots = a
    if index == 1:
        numplots = b
    elif index == 2:
        numplots = c
        
    for i in range(numplots):
        if index == 0:
            data = recon[i,:,:]
        elif index == 1:
            data = recon[:,i,:]
        else:
            data = recon[:,:,i]
        plt.figure()
        plt.imshow(data, interpolation = interpolation)
        plt.axis('off')
        if savefig:
            plt.savefig(os.path.join(outdir, '%s_axis%d_%d.png'%(smp,index,i)),
                        dpi = 300, transparent = True)
    
        
        
def plotReconSlice(recon1, recon2, slice_dir, interpolation = None, 
                   savefig = None, outDir = None):
    num_slice = recon1.shape[slice_dir]
    
    for i in range(num_slice):
        if slice_dir == 0:
            data1 = recon1[i,:,:]
            data2 = recon2[i,:,:]
        elif slice_dir == 1:
            data1 = recon1[:,i,:]
            data2 = recon2[:,i,:]
        else:
            data1 = recon1[:,:,i]
            data2 = recon2[:,:,i]
            
        f = plt.figure(figsize=(3,3))
        a = plt.imshow(data1, cmap = plt.cm.get_cmap('binary_r'), alpha = 1, 
                        interpolation = interpolation)
        b = plt.imshow(data2, cmap = plt.cm.get_cmap('Blues'), alpha = 0.5, 
                        interpolation = interpolation)
#        a = plt.imshow(data1, cmap = plt.cm.get_cmap('binary_r'), alpha = 1, 
#                       vmin = 0.03, vmax = 0.1, interpolation = interpolation)
#        b = plt.imshow(data2, cmap = plt.cm.get_cmap('Blues'), alpha = 0.5, 
#                       vmin = 0.1, vmax = 0.5, interpolation = interpolation)
        plt.axis('off')
        ax = f.gca()
        
        cax1 = f.add_axes([1, 0.05, 0.03, 0.9])
        cax2 = f.add_axes([1.2, 0.05, 0.03, 0.9])
        cbar = f.colorbar(a, cax=cax1, orientation='vertical', ticks=[0.03, 0.10])
        cbar2 = f.colorbar(b, cax=cax2, orientation='vertical', ticks=[0.10, 0.50])
        cbar.set_label('S', color = plt.cm.get_cmap('binary_r')(0.1), labelpad = -20)
        cbar2.set_label('Cd', color=plt.cm.get_cmap('Blues')(0.9),  labelpad = -10)
        plt.tight_layout()
                
        if savefig:
            plt.savefig(os.path.join(outDir, 'slice%d_%s_%d.png'%(slice_dir, interpolation, i)),
                        dpi = 300, bbox_inches = "tight")
            
def framesToVideo(image_folder, smp_name, splitIndex, output_dir = None, 
                  algorithm = None, n_proj = None):
    video_name = '%s_%s_%s.mov'%(algorithm, n_proj, smp_name)
    images = [img for img in os.listdir(image_folder) if smp_name in img]
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape
    fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v') # note the lower case
    
    if output_dir is None:
        output_dir = image_folder
    
    video = cv2.VideoWriter(os.path.join(output_dir, video_name), fourcc, 5, (width,height))
    images.sort(key = lambda a: int(a.split('_')[splitIndex][:-4]))
    
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))
    
    cv2.destroyAllWindows()
    video.release()