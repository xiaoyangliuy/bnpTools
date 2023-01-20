#!/usr/bin/env python
# coding: utf-8

import astra
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import os, sys, h5py, json
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import tifffile

plt.gray()
get_ipython().run_line_magic('matplotlib', 'inline')


# Set up astra to reconstruct a cross-section slice from a sinogram
def sinoRecon(sinogram, angles, w_pixel, algorithm, n_iter = None, ycenter = 0, xcenter = 0, volscale = 1):
    n_angles, n_y = sinogram.shape
    assert (n_angles == len(angles)), 'Dim of sinogram and angle does not match. Sinogram should have shape of (angles, voxels)'
    
    proj_geom = astra.create_proj_geom('parallel', w_pixel, n_y, angles)
    proj_geom = astra.geom_postalignment(proj_geom, [ycenter, xcenter])
        
    vol_geom = astra.create_vol_geom(int(n_y*volscale), int(n_y*volscale))
    proj_id = astra.create_projector('linear', proj_geom, vol_geom)
    
    sinogram_id = astra.data2d.link('-sino', proj_geom, sinogram)
    recon_id = astra.data2d.create('-vol', vol_geom)
    
    cfg = astra.astra_dict(algorithm)
    cfg['ProjectorId'] = proj_id
    cfg['ProjectionDataId'] = sinogram_id
    cfg['ReconstructionDataId'] = recon_id
    
    if n_iter is None:
        n_iter = 1
    algorithm_id = astra.algorithm.create(cfg)
    astra.algorithm.run(algorithm_id, n_iter)
    return astra.data2d.get(recon_id)


# In[ ]:


def orderProjAxis(proj, angle_axis, col_axis, row_axis, elm_axis):
    ndim = len(proj.shape)
    newProj = proj.copy()
    
    if elm_axis is None:
        if ndim > 3:
            raise ValueError('Need to specify axis dimension of elemental channel')
        elif ndim == 3:
            newProj = np.expand_dims(proj, axis=0)
            elm_axis = 0
    else:
        if elm_axis > ndim:
            raise ValueError('Elm axis is larger than array dimension')

    # move elm-, angle-, col- and row- axis to the designated position
    org = {'elm':elm_axis, 'angle':angle_axis, 'col':col_axis, 'row': row_axis}
    org_sort = {k: v for k, v in sorted(org.items(), key=lambda item: item[1], reverse = True)}
    goal = {'elm':0, 'col':1, 'angle':2, 'row': 3}
    for i, k in enumerate(org_sort.keys()):
        newProj = np.moveaxis(newProj, org_sort[k], goal[k])
        dict_k = list(org_sort.keys())[i:]
        if org_sort[k] != goal[k]:
            for j in dict_k:
                    org_sort[j] += 1
    return newProj


# Iterate through all sinogram slice and return the reconstructed volume 
def recon(proj, angles, angle_axis, col_axis, row_axis, elm_axis, w_pixel, algorithm, 
          n_iter = None, ycenter = 0, xcenter = 0, volscale = 1):
    oProj = orderProjAxis(proj, angle_axis = angle_axis, col_axis = col_axis, 
                          row_axis=row_axis, elm_axis = elm_axis)
    n_elm, n_sino, n_angle, n_row = oProj.shape
    r = np.zeros((n_elm, n_sino, int(n_row * volscale), int(n_row * volscale)))
    for i in range(n_elm):
        for j in range(n_sino):
            sino = np.array(oProj[i,j,:,:], order = 'C', dtype='float32')
            vslice = sinoRecon(sino, angles, w_pixel, algorithm, n_iter = n_iter, 
                               ycenter = ycenter, xcenter = xcenter, volscale = volscale)
            r[i,j,:,:] = vslice
    return {'proj_input':oProj, 'recon':r}

def sinoRecon_2DCUDA(sinogram, angles, w_pixel, algorithm='SIRT_CUDA', n_iter = 1, 
                     ycenter = 0, xcenter = 0, volscale = 1, miniConstraint = 0):
    n_angles, n_y = sinogram.shape
    assert (n_angles == len(angles)), 'Dim of sinogram and angle does not match. Sinogram should have shape of (angles, voxels)'
    
    proj_geom = astra.create_proj_geom('parallel', w_pixel, n_y, angles)
    proj_geom = astra.geom_postalignment(proj_geom, [ycenter, xcenter])
    
    vol_geom = astra.create_vol_geom(int(n_y*volscale), int(n_y*volscale))
    proj_id = astra.create_projector('cuda', proj_geom, vol_geom)
    
    sinogram_id = astra.data2d.link('-sino', proj_geom, sinogram)
    recon_id = astra.data2d.create('-vol', vol_geom)
    
    cfg = astra.astra_dict(algorithm)
    cfg['ProjectorId'] = proj_id
    cfg['ProjectionDataId'] = sinogram_id
    cfg['ReconstructionDataId'] = recon_id
    cfg['MinConstraint'] = miniConstraint
    
    algorithm_id = astra.algorithm.create(cfg)
    
    if n_iter is None:
        astra.algorithm.run(algorithm_id)
    else:
        astra.algorithm.run(algorithm_id, n_iter)
    return astra.data2d.get(recon_id)


# Iterate through all sinogram slice and return the reconstructed volume 
def recon_CUDA(proj, angles, algorithm, angle_axis=2, col_axis = 1, row_axis = 3, elm_axis = 0, w_pixel = 1,  
          n_iter = 300, ycenter = 0, xcenter = 0, volscale = 1):
    oProj = orderProjAxis(proj, angle_axis = angle_axis, col_axis = col_axis, 
                          row_axis=row_axis, elm_axis = elm_axis)
    n_elm, n_sino, n_angle, n_row = oProj.shape
    r = np.zeros((n_elm, n_sino, int(n_row * volscale), int(n_row * volscale)))
    for i in range(n_elm):
        for j in range(n_sino):
            sino = np.array(oProj[i,j,:,:], order = 'C', dtype='float32')
            vslice = sinoRecon_2DCUDA(sino, angles, w_pixel, algorithm, n_iter = n_iter, 
                               ycenter = ycenter, xcenter = xcenter, volscale = volscale)
            vslice[vslice < 0] = 0
            r[i,j,:,:] = vslice
    return {'proj_input':oProj, 'recon':r}


#%% Select limited projections
def numProjIndex(thetas, numproj):
    projIdx = None
    if numproj <= len(thetas):
        projIdx = np.arange(0, 180, 180/numproj, dtype='int')
    else:
        raise ValueError('numproj is larger than theta size')
    return projIdx


# In[ ]:


# %% Plot the reproj images 
def plotRepro(numprojs, selReconData, proj, proj_axis = None, sino = None, 
              sino_slice = None, cmax=None, coltitle=None, figsize = None):
    labels = ['P', 'Zn', 'Fe']
    cmap = ['inferno', 'viridis', 'cividis']
    tpos = [-30, -40, -40]
    
    if proj_axis is None:
        proj_axis = 2
        
    if proj_axis != 0:
#         figsize = (8,4)
        fig, axes = plt.subplots(3,len(numprojs)+1, figsize=figsize)
        cbar_s = 0.98
#         cmax = [1.5, 0.05, 0.0055]
        for i in range(3):
            data = proj[0,:,:,i]
            m = axes[i,-1].imshow(data, cmap[i], vmin = 0, vmax = cmax[i])
            axes[i,-1].set_xticks([])
            axes[i,-1].set_yticks([])
            cbar = fig.colorbar(m, ax = axes[i,-1],shrink=cbar_s)
        axes[0,-1].set_title('Proj')
        
    elif (proj_axis == 0) | (sino is True):
#         figsize = (9.5,4)
        fig, axes = plt.subplots(3,len(numprojs), figsize=figsize)
        cbar_s = 0.98
#         cmax = [2.2, 0.12, 0.0055]  
        
    for i, k in enumerate(selReconData.keys()):
        for j, e in enumerate(labels):
            if sino:
                data = selReconData[k]['proj_input'][j,sino_slice[j],:,:]
                # axes[j,i].set_aspect('auto')
                vmax = cmax[j]
            else:
                data = np.sum(selReconData[k]['recon'][j,:,:,:], axis=proj_axis).T
                vmax = cmax[j]
                
            m = axes[j,i].imshow(data, cmap[j], vmin = 0, vmax=vmax)
            axes[j,i].axis('off')
            axes[j,i].set_xticks([])
            axes[j,i].set_yticks([])
#             cbar = fig.colorbar(m, ax = axes[j,i], shrink=cbar_s)
            if i == 0:
                axes[j,i].text(tpos[j],70,labels[j])
            if (j == 0) & (coltitle is None) :
                axes[j,i].set_title(str(selReconData[k]['proj_input'].shape[2]))
            elif j == 0:
                axes[j,i].set_title(coltitle[i])
            if (proj_axis == 0) & (i == len(numprojs)-1):
                cbar = fig.colorbar(m, ax = axes[j,i],shrink=cbar_s)
    plt.tight_layout()
    return fig, axes


# In[420]:


# Plot reproj lineprofile
def plotReproLine(numprojs, selReconData, proj, proj_axis=None, ymax = None):
    labels = ['P', 'Zn', 'Fe']
    cmap = ['Reds', 'Greens', 'Blues']
    label_c = ['r', 'g', 'b']
    xlim = [85, 85, 75]
    ymin = [0.6, 0.02, -2e-6]
#     ymin = [0.75, 0.025, -1e-6]
    
    
    
    fig, axes = plt.subplots(len(labels), 1, figsize = (5,7), sharex=True)
    
    if proj_axis is None:
        proj_axis = 2
    
    for j, e in enumerate(labels):
        cidx = np.linspace(0.9,0.2,len(selReconData))
        axins2 = inset_axes(axes[j], width="30%", height="40%")
        for i, k in enumerate(selReconData.keys()):
            p = np.sum(selReconData[k]['recon'][j,:,:,:], axis=proj_axis).T
            lineprofile = np.sum(p, axis = 1) / p.shape[1]
            x = np.arange(0,len(lineprofile))
            axes[j].plot(x, lineprofile,
                         color = plt.get_cmap(cmap[j])(cidx[i]), label=k)
            axins2.plot(x, lineprofile, color = plt.get_cmap(cmap[j])(cidx[i]))
            
            
        p_proj = proj[0,:,:,j]
        pline = np.sum(p_proj, axis = 1) / p_proj.shape[1]
        axes[j].plot(np.arange(0,len(pline)), pline,
                         color = 'k', label='Proj', alpha = 0.8)
        axins2.plot(np.arange(0,len(pline)), pline,
                         color = 'k', label='Proj', alpha = 0.8)
        y_min, y_max = axes[j].get_ylim()
        axins2.set_ylim((ymin[j], y_max))
        axins2.set_yticks([])
        
        if ymax is not None:
            y_max = ymax[j]
            axes[j].set_ylim((y_min, y_max))
            
        axes[j].text(0, y_max*0.85, e, color = label_c[j])
        axes[j].set_ylabel('Average Pixel Intensity (a.u.)')
        axins2.set_xlim((50, xlim[j]))
        axes[j].legend(frameon=False, bbox_to_anchor=(1, 0.3, 0.5, 0.5))
    
    axes[j].set_xlabel('Position (pixel)')
    return fig



def misalign_proj(proj, misalign_ratio, axis = None):
    
    # proj is a np array with dimension of [angle, row, col, elm]
    # axis is the misaligned axis. axis = 1, x-axis; axis = 2, y-axis, axis = 3, x- and y- axes
    n_angle = proj.shape[0]
    n_proj = int(n_angle * misalign_ratio)
    
    # select random projection and shifts
    seed = 40
    np.random.seed(seed)
    nproj = np.random.randint(n_angle, size=(n_proj))
    
    if axis is None:
        axis = 1
    
    if axis <= 2:
        n_pix = proj.shape[axis]
        shifts = np.random.randint(-n_pix * 0.2, high = n_pix * 0.2, size=(n_proj,1))
        s_axis = axis - 1
    elif axis == 3:
        pix1, pix2 = proj[0,:,:,0].shape
        s1 = np.random.randint(-pix1 * 0.2, high = pix1 * 0.2, size=(n_proj))
        s2 = np.random.randint(-pix2 * 0.2, high = pix2 * 0.2, size=(n_proj))
        shifts = np.vstack((s1,s2)).T
        s_axis = (0, 1)
    else:
        raise ValueError('Value of n_axis is too large. Expecting n_axis (number of misaligned axes to be less than or equal to 3)')

    nonalign_proj = proj.copy()
    for n_, s_ in zip(nproj, shifts):
        temp = nonalign_proj[n_,:,:,:].copy()
        nonalign_proj[n_,:,:,:] = np.roll(temp, tuple(s_), axis = s_axis)
        
    return nonalign_proj


def write_tiff(data, fname='tmp/data', digit=None, ext='tiff'):
    """
    Write image data to a tiff file.

    Overwrite existing data and infer data-type from the data.

    Parameters
    ----------
    data : ndarray
        Array data to be saved.
    fname : str
        File name to which the data is saved. ``.tiff`` extension
        will be appended if it doesn't already have one.
    digit : int
        Append this number to fname using a folder e.g. {fname}/{digit}.{ext}
    """
    # Add the extension and digit.
    if digit is not None:
        fname = os.path.join(fname, str(digit))
    if not str(fname).endswith(ext):
        fname = ".".join([fname, ext])
    # Convert to absolute path.
    fname = os.path.abspath(fname)
    # Create the directory if it doesn't exist.
    dname = os.path.dirname(os.path.abspath(fname))
    if not os.path.exists(dname):
        os.makedirs(dname)
    # Save the file.
    tifffile.imsave(fname, data)
    
    
#%% Export data in .h5 format for visualization in Tomviz
def makeTomvizH5(fname, data, scans, theta, xval, yval):
    dt = h5py.string_dtype(encoding='utf-8')
    with h5py.File(fname, 'w') as f:
        g1 = f.create_group('exchange')
        g1.create_dataset('data',data = data)
        g1.create_dataset('scanName', data = scans)
        g1.create_dataset('theta', data = theta)
        g1.create_dataset('x_axis', data = xval)
        g1.create_dataset('y_axis', data = yval)
        
def loadAlignedDataFromJson(fpath):
    with open(fpath) as json_file:
        data = json.load(json_file)

    elms = list(data.keys())[3:]
    angles = np.array(data['angles'], dtype = 'float')
    a, h, w = np.array(data[elms[0]]).shape
    elmdata = np.zeros((len(elms), a, h, w))
    
    for i, e in enumerate(elms):
        elmdata[i,...] = np.array(data[e])
    elmdata = np.moveaxis(elmdata, 1, 2)
    p_angles = angles * np.pi / 180
    
    return {'elms':elms, 'elmdata':elmdata, 'angles_deg':angles, 'angle_rad':p_angles, 'stepsize':data['stepsize']}





