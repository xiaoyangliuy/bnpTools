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



# Show the overview map to reflect scan coordinate
def plotoverview(df, sc, ax = None):
    if ax is None:
        fig = plt.figure(figsize = (3, 3))
        ax = plt.gca()
        
    box_params = []
    for i, row in df.iterrows():
        x = row['xval']
        y = row['yval']
        corner = (min(x), min(y))
        x_width = max(x) - min(x)
        y_width = max(y) - min(y)
        scan_area = (x_width * y_width)
        scan = row['scan_num']
        box_params.append([scan_area, corner, x_width, y_width, scan])

    for i, box in enumerate(box_params):
        if sc == box[-1]:
            alpha = 0.8
        else:
            alpha = 0.15
        color = plt.cm.get_cmap('Pastel2')(i)
        hr = Rectangle(box[1], box[2], box[3], picker = True, facecolor = color, 
                       edgecolor = [0, 0, 0], label = 'scan:%s'%box[-1], alpha = alpha)
        ax.add_patch(hr)

    ax.autoscale(enable = True)
    
    return ax



# K-mean cluster analysis
def kmean_analysis(n_clusters, data, random_state, sigma = None, cval = None,
                   plotoption = None, savefig = None, fname = None):
    
    data[np.isnan(data)] = 1e-5
    data[np.isinf(data)] = 1e-5
    
    if sigma is None: 
        data_blur = data
    else:
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


def singleChArray(df, ch_name):
    ysize, xsize = df[ch_name][0].shape
    data = np.zeros([1, len(df), ysize, xsize])
    for i, row in df.iterrows():
        d = row[ch_name]
        d[np.isnan(d)] = 1e-3
        d[np.isinf(d)] = 1e-3
        try:
            data[0,i,:,:] = d
        except:
            print(d.shape, i, ysize, xsize)
    return data

def multiChArray(data, elms):
    ss = np.asarray([m_.shape for m_ in data[elms[0]]])
    ysize = np.max(ss[:, 0])
    xsize = np.max(ss[:, 1])
    
    elmdata = np.zeros([len(elms), len(data), ysize, xsize])
    for i, elm in enumerate(elms):
        n = 0
        for j, row in data.iterrows():
            d = row[elm][:]
            d[np.isnan(d) | np.isinf(d)] = 1e-3
            if d.size == (ysize * xsize):
                elmdata[i,n,:,:] = d
            else:
                dy = (ysize - d.shape[0])//2
                dx = (xsize - d.shape[1])//2
                elmdata[i, n, dy:d.shape[0]+dy, dx:d.shape[1]+dx] = d
            n += 1
    return elmdata


# def multiChArray(df, ch_names, dfIndex = None):
#     ysize, xsize = df[ch_names[0]][0].shape
    
#     if not dfIndex:
#         data = np.zeros([len(ch_names), len(df), ysize, xsize])
#         for i, elm in enumerate(ch_names):
#             n = 0
#             for j, row in df.iterrows():
#                 d = row[elm][:].copy()
#                 d[np.isnan(d)] = 1e-3
#                 d[np.isinf(d)] = 1e-3
#                 data[i,n,:,:] = d
#                 n += 1
#     else:
#         data = np.zeros([len(ch_names), len(dfIndex), ysize, xsize])
#         for i, elm in enumerate(ch_names):
#             n = 0
#             for idx, j in enumerate(dfIndex):
#                 d = df.iloc[j][elm].copy()
#                 d[np.isnan(d)] = 1e-3
#                 d[np.isinf(d)] = 1e-3
#                 data[i,n,:,:] = d
#                 n += 1
#     return data

def getSinoData(projData,axisIndex):
    # assume projData is 3D array with dimension of [theta,y,x]
    sino = np.zeros((projData.shape[1],projData.shape[0]))
    for i in range(len(projData)):
        sino[:,i] = -np.log(np.sum(projData[i,:,:],axis = axisIndex))
    return sino


#%% Assess quality of reconstruction
def assessRecon(recon, data, thetas, mid_indx, show_plots=False):
    zero_index = np.argmin(abs(thetas))
    num_slices = recon.shape[0]
    width = recon.shape[1]
    reprojection = np.zeros([num_slices, width])
    tmp = np.zeros([num_slices, width])
   
    # get recon reporjection for slice i and take the difference with data projection (at angle ~=0).
    for i in range(num_slices):
        reprojection[i] = np.sum(recon[i], axis=0)
        if data[zero_index, i].max() == 0:
            tmp[i] = np.zeros(width)
        else:
            #projection for 0 angle at row i / projection for angle 0 at row i, 
            #mximum value / maximum value of recon for row i
            #projection row / (proj max / reproj max)
            tmp[i] = data[zero_index, i] / (data[zero_index, i].max() / np.sum(recon[i], axis=0).max())
    #normalizing projectios against reconstruction side length, so plot appears withing image.
    projection = tmp*width/tmp.max()
    #normalizing reprojectios against reconstruction side length, so plot appears withing image.
    nonNormrepro = reprojection
    reprojection = reprojection*width/reprojection.max()
  
    projection_xSection = tmp[mid_indx]*width/tmp[mid_indx].max()
    reprojection_xSection = reprojection[mid_indx]*width/reprojection[mid_indx].max()
    #difference between reporjection and original projection at angle == 0
    # err = tmp - reprojection/reprojection
    err = projection - reprojection
    #mean squared error
    mse = (np.square(err)).mean(axis=None)
    if show_plots:
        plt.figure()
        plt.imshow(recon[mid_indx], origin='lower')
        plt.plot(projection_xSection)
        plt.plot(reprojection_xSection)
        plt.legend(('projection', 'reprojection'), loc=1)
        plt.title("MSE:{}".format(np.round(mse, 4)))
        plt.figure()
        plt.imshow(projection)
        plt.title("projection")
        plt.colorbar()
        plt.figure()
        plt.imshow(reprojection)
        plt.title("reprojection")
        plt.colorbar()
        plt.show()
        
    return err, mse, nonNormrepro

#%% Export data in .h5 format for visualization in Tomviz
def makeTomvizH5(fname, data, scans, theta, xval, yval):
    dt = h5py.string_dtype(encoding='utf-8')
    with h5py.File(fname, 'w') as f:
        g1 = f.create_group('exchange')
        g1.create_dataset('data',data = data)
        g1.create_dataset('scanName', data = scans)
        g1.create_dataset('theta', data = theta)
        g1.create_dataset('x_axis', data = xval, dtype = dt)
        g1.create_dataset('y_axis', data = yval, dtype = dt)
        
        
#%% Feature matching
def siftFeatureDetect(im1,im2, sigma = None, nOct = None, im1mask=None, 
                      im2mask = None, limit = None, axes = None, plotFig = True):
    sift = cv2.SIFT_create(sigma=sigma, nOctaveLayers = nOct)
    kp1, des1 = sift.detectAndCompute(im1,im1mask)
    kp2, des2 = sift.detectAndCompute(im2,im2mask)
    
    FLANN_INDEX_KDTREE = 0
    index_params = dict(algorithm = FLANN_INDEX_KDTREE, trees = 5)
    search_params = dict(checks = 500)
    
    flann = cv2.FlannBasedMatcher(index_params, search_params)
    matches = flann.knnMatch(des1,des2,k=2)
    im2nolines = im2.copy()
        
    # store all the good matches as per Lowe's ratio test.
    good = []
    if limit is None:
        limit = 0.7
    for m,n in matches:
        if m.distance < limit*n.distance:
            good.append(m)
    
    MIN_MATCH_COUNT = 4
    dst = None
    M = None
    if len(good)>=MIN_MATCH_COUNT:
        src_pts = np.float32([ kp1[m.queryIdx].pt for m in good ]).reshape(-1,1,2)
        dst_pts = np.float32([ kp2[m.trainIdx].pt for m in good ]).reshape(-1,1,2)
    
        M, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 5.0, 
                                     maxIters = 5000, confidence=0.99)
        matchesMask = mask.ravel().tolist()
    
        h,w = im1.shape
        pts = np.float32([ [0,0],[0,h-1],[w-1,h-1],[w-1,0] ]).reshape(-1,1,2)
        dst = cv2.perspectiveTransform(pts,M)
    
        im2 = cv2.polylines(im2,[np.int32(dst)],True,255,3, cv2.LINE_AA)
    
    else:
        print ("Not enough matches are found - %d/%d" % (len(good),MIN_MATCH_COUNT))
        matchesMask = None
        
    draw_params = dict(matchColor = (0,255,0), # draw matches in green color
                        singlePointColor = None,
                        matchesMask = matchesMask, # draw only inliers
                        flags = 2)
    
    img3 = cv2.drawMatches(im1,kp1,im2,kp2,good,None,**draw_params)
    # img3 = cv2.drawMatchesKnn(im1,kp1,im2,kp2,good,outImg=None,flags=2)
    img2 = cv2.drawKeypoints(im2nolines, kp2, None, color = (0,255,0))
    img1 = cv2.drawKeypoints(im1.copy(), kp1, None, color = (0,255,0))
    
    if axes is not None:
        axes.imshow(img3, 'gray')
        axes.axis('off')
    if plotFig:
        plt.imshow(img3,'gray')
        plt.axis('off')

    return {'combined':img3, 'kpts_img2':img2, 'kpts_img1':img1, 'M':M, 'dstPts':dst,
           'im2_box':im2, 'good':good}


#%% Feature matching (limited to rigid 2D transform)
def siftFeatureDetect_AffineAartial2D(im1,im2, sigma = None, nOct = None, im1mask=None, 
                      im2mask = None, limit = None, axes = None, plotFig = True):
    sift = cv2.SIFT_create(sigma=sigma, nOctaveLayers = nOct)
    kp1, des1 = sift.detectAndCompute(im1,im1mask)
    kp2, des2 = sift.detectAndCompute(im2,im2mask)
    
    FLANN_INDEX_KDTREE = 0
    index_params = dict(algorithm = FLANN_INDEX_KDTREE, trees = 5)
    search_params = dict(checks = 500)
    
    flann = cv2.FlannBasedMatcher(index_params, search_params)
    matches = flann.knnMatch(des1,des2,k=2)
    im2nolines = im2.copy()
        
    # store all the good matches as per Lowe's ratio test.
    good = []
    if limit is None:
        limit = 0.7
    for m,n in matches:
        if m.distance < limit*n.distance:
            good.append(m)
    
    MIN_MATCH_COUNT = 3
    dst = None
    M = None
    if len(good)>=MIN_MATCH_COUNT:
        src_pts = np.float32([ kp1[m.queryIdx].pt for m in good ]).reshape(-1,1,2)
        dst_pts = np.float32([ kp2[m.trainIdx].pt for m in good ]).reshape(-1,1,2)
    
        M, mask = cv2.estimateAffinePartial2D(src_pts, dst_pts, M, int(cv2.RANSAC), 5.0, 
                                     maxIters = 5000, confidence=0.99)
        matchesMask = mask.ravel().tolist()
    
        h,w = im1.shape
        pts = np.float32([ [0,0],[0,h-1],[w-1,h-1],[w-1,0] ]).reshape(-1,1,2)
        dst = cv2.transform(pts, M)
#         dst = cv2.affineTransform(pts, M)
#         dst = cv2.perspectiveTransform(pts,M)
    
        im2 = cv2.polylines(im2,[np.int32(dst)],True,255,3, cv2.LINE_AA)
    
    else:
        print ("Not enough matches are found - %d/%d" % (len(good),MIN_MATCH_COUNT))
        matchesMask = None
        
    draw_params = dict(matchColor = (0,255,0), # draw matches in green color
                        singlePointColor = None,
                        matchesMask = matchesMask, # draw only inliers
                        flags = 2)
    
    img3 = cv2.drawMatches(im1,kp1,im2,kp2,good,None,**draw_params)
    # img3 = cv2.drawMatchesKnn(im1,kp1,im2,kp2,good,outImg=None,flags=2)
    img2 = cv2.drawKeypoints(im2nolines, kp2, None, color = (0,255,0))
    img1 = cv2.drawKeypoints(im1.copy(), kp1, None, color = (0,255,0))
    
    if axes is not None:
        axes.imshow(img3, 'gray')
        axes.axis('off')
    if plotFig:
        plt.imshow(img3,'gray')
        plt.axis('off')

    return {'combined':img3, 'kpts_img2':img2, 'kpts_img1':img1, 'M':M, 'dstPts':dst,
           'im2_box':im2}

def getNorm8U(img):
    imgNorm = 255*((img-np.amin(img))
                  /(np.amax(img)-np.amin(img)))
    imgNorm = imgNorm.astype('uint8')
    return imgNorm

#%% Order pts from perspective transform
def order_points(pts):
    rect = np.zeros((4, 2), dtype = "float32")
    s = pts.sum(axis = 1)
    rect[0] = pts[np.argmin(s)]
    rect[2] = pts[np.argmax(s)]
    diff = np.diff(pts, axis = 1)
    rect[1] = pts[np.argmin(diff)]
    rect[3] = pts[np.argmax(diff)]
    return rect

def four_point_transform(image, pts, M=None):
#     rect = order_points(pts)
    rect = order_points(pts)
    (tl, tr, br, bl) = rect
    # compute the width of the new image, which will be the
    # maximum distance between bottom-right and bottom-left
    # x-coordiates or the top-right and top-left x-coordinates
    widthA = np.sqrt(((br[0] - bl[0]) ** 2) + ((br[1] - bl[1]) ** 2))
    widthB = np.sqrt(((tr[0] - tl[0]) ** 2) + ((tr[1] - tl[1]) ** 2))
    maxWidth = max(int(widthA), int(widthB))
    # compute the height of the new image, which will be the
    # maximum distance between the top-right and bottom-right
    # y-coordinates or the top-left and bottom-left y-coordinates
    heightA = np.sqrt(((tr[0] - br[0]) ** 2) + ((tr[1] - br[1]) ** 2))
    heightB = np.sqrt(((tl[0] - bl[0]) ** 2) + ((tl[1] - bl[1]) ** 2))
    maxHeight = max(int(heightA), int(heightB))
    # now that we have the dimensions of the new image, construct
    # the set of destination points to obtain a "birds eye view",
    # (i.e. top-down view) of the image, again specifying points
    # in the top-left, top-right, bottom-right, and bottom-left
    # order
    dst = np.array([
        [0, 0],
        [maxWidth - 1, 0],
        [maxWidth - 1, maxHeight - 1],
        [0, maxHeight - 1]], dtype = "float32")
    # compute the perspective transform matrix and then apply it
    M = cv2.getPerspectiveTransform(rect, dst)
    warped = cv2.warpPerspective(image, M, (maxWidth, maxHeight))
    # return the warped image
    return warped, M


def fft_radical_average(rcn):
    f = np.fft.fft2(rcn)
    fshift = np.fft.fftshift(f)
    magnitude = np.abs(fshift)
    rows, cols = rcn.shape
    midrow = rows/2 + 1
    midcol = cols/2 + 1
    maxr = np.ceil(np.sqrt((rows)**2 + (cols)**2))
    radical = np.zeros(int(maxr))
    count = np.zeros(int(maxr))
    for c in range(cols):
        for r in range(rows):
            radius = np.sqrt((r - midrow)**2 + (c - midcol)**2)
            thisIdx = int(np.ceil(radius) + 1)
            radical[thisIdx] = radical[thisIdx] + magnitude[r, c]
            count[thisIdx] = count[thisIdx] + 1

    radicalProfile = radical / count
    return radicalProfile

def azimuthalAverage(image, center=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr

    return radial_prof
