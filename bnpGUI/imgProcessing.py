#!/home/beams/USERBNP/.conda/envs/py36/bin/python
import numpy as np
import h5py, os
from sklearn.cluster import KMeans
from skimage.measure import regionprops
from skimage.draw import rectangle, polygon_perimeter
from skimage.transform import resize 
from scipy.fftpack import fft2, ifft2
import skimage.filters
import matplotlib.patches as mpatches
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage import io
import cv2 as cv
from scipy.ndimage import shift, center_of_mass
from scipy import ndimage
import epics as PV
import matplotlib.patches as patches

def estShifts(refSc, currSc, fpath, elm):
    ref, ref_x_pos, ref_y_pos = getElmMap(os.path.join(fpath, refSc), elm)
    mov, mov_x_pos, mov_y_pos = getElmMap(os.path.join(fpath, currSc), elm)
    ref_resize = resize(ref, mov.shape)
    s1, s2 = phaseCorrelate(ref_resize, mov)
    fig = plt.figure()
    plt.imshow(ref_resize, cmap='Greens', alpha = 0.5)
    plt.imshow(np.roll(mov.copy(), (int(s1), int(s2))), cmap='Reds', alpha = 0.5)
    ss_x = np.diff(mov_x_pos)[0]*s1
    ss_y = np.diff(mov_y_pos)[0]*s2
    plt.title('x_shift_um:%.2f, y_shift_um:%.2f, x_shift_pixel:%d, y_shift_pixel:%d'%(ss_x, ss_y, s1, s2))
    return ss_x, ss_y, fig

def getROIcoordinate(fname, elm, savefig = True, figpath = None, n_cluster = 2, sel_cluster = 1):
    elmmap, x_pos, y_pos = getElmMap(fname, elm)
    kmeanMap = kmean_analysis(n_cluster,elmmap, 42, plotoption = True)
    region_prop = regionprops(np.array(kmeanMap[0]==sel_cluster, dtype='int'))
    region_bbox = region_prop[0].bbox
    fig = plotBBox(elmmap, region_bbox, x_pos, y_pos)
        
    width = x_pos[region_bbox[3]-1]- x_pos[region_bbox[1]]
    height = y_pos[region_bbox[2]-1]- y_pos[region_bbox[0]]
    
    new_x = (x_pos[region_bbox[3]-1]- x_pos[region_bbox[1]])/2 + x_pos[region_bbox[1]]
    new_y = (y_pos[region_bbox[2]-1]- y_pos[region_bbox[0]])/2 + y_pos[region_bbox[0]]
    
    if savefig:
        fig.savefig(figpath, dpi=100, transparent=True)
        
    return new_x, new_y, width, height

def getROIcoordinate_data(elmmap, x_pos, y_pos, savefig = True, figpath = None, n_cluster = 2, sel_cluster = 1):
    kmeanMap = kmean_analysis(n_cluster,elmmap, 42, plotoption = False)
    region_prop = regionprops(np.array(kmeanMap[0]==sel_cluster, dtype='int'))
    region_bbox = region_prop[0].bbox
    fig = plotBBox(elmmap, region_bbox, x_pos, y_pos)
        
    width = x_pos[region_bbox[3]-1]- x_pos[region_bbox[1]]
    height = y_pos[region_bbox[2]-1]- y_pos[region_bbox[0]]
    
    new_x = (x_pos[region_bbox[3]-1]- x_pos[region_bbox[1]])/2 + x_pos[region_bbox[1]]
    new_y = (y_pos[region_bbox[2]-1]- y_pos[region_bbox[0]])/2 + y_pos[region_bbox[0]]
    
    if savefig:
        fig.savefig(figpath, dpi=100, transparent=True)
        
    return new_x, new_y, width, height

def checkROIIntensity(fname, elm, coordinates):
    elmmap, x_pos, y_pos = getElmMap(fname, elm)
    
    
def plotBBox(elmmap, box, x_pos, y_pos):
    minr, minc, maxr, maxc = box
    x_st = x_pos[box[1]]
    y_st = y_pos[box[0]]
    w = y_pos[box[2]-1] - y_pos[box[0]]
    h = x_pos[box[3]-1] - x_pos[box[1]]
    rect = mpatches.Rectangle((x_st, y_st), h, w,
            fill=False, edgecolor='red', linewidth=2)
    fig = plt.figure()
    ax = fig.gca()
    ax.pcolor(x_pos, y_pos, elmmap, cmap='gray', shading='auto')
    ax.add_patch(rect)
    return fig
'''
def getElmMap(fname, elm):
    elmmap = []
    with h5py.File(fname, 'r') as dat:
        xrfdata = dat['/MAPS/XRF_roi'][:]
        try:
            ch_idx = dat['/MAPS/channel_names'][:].astype('U13').tolist().index(elm)
            elmmap = xrfdata[ch_idx,:,:]
            x_pos = dat['/MAPS/x_axis'][:]
            y_pos = dat['/MAPS/y_axis'][:]
        except:
            raise ValueError('Invalid element! %s is not in channel list'%(elm))
    return elmmap, x_pos, y_pos
'''
#---------------------------xyl----------------------------
def getElmMap(elm,flag, test_folder=None, test_sc_range=None,fname=''):
    #elmmap = []
    if flag == 0:
        with h5py.File(fname, 'r') as dat:
            xrfdata = dat['/MAPS/XRF_roi'][:]
            try:
                channel_names = np.char.decode(dat['MAPS']['channel_names'][:])
                elm_idx = np.where(channel_names==elm)[0][0]
                suc = True    
                elmmap = xrfdata[elm_idx,:,:]
                x_pos = dat['/MAPS/x_axis'][:]
                y_pos = dat['/MAPS/y_axis'][:]
                theta = round(float(dat['MAPS/extra_pvs'][1,8].decode('utf-8')))  #get angle int info
            except IndexError:
                print('Invalid element! %s is not in channel list'%(elm))
                suc = False
                elmmap = None
                x_pos = 0        
                y_pos = 0      
                theta = round(float(dat['MAPS/extra_pvs'][1,8].decode('utf-8')))  #get angle int info
        f = 0
    else:
        scan_num = test_sc_range[getElmMap.counter]
        scan_num_formatted = '{:04d}'.format(scan_num)
        f = f'{test_folder}bnp_fly{scan_num_formatted}.mda.h5'
        print(f'test_image:{f}')
        with h5py.File(f, 'r') as dat:
            xrfdata = dat['/MAPS/XRF_roi'][:]
            channel_names = np.char.decode(dat['MAPS']['channel_names'][:])
            elm_idx = np.where(channel_names==elm)[0][0]
            elmmap = xrfdata[elm_idx,:,:]
            x_pos = dat['/MAPS/x_axis'][:]
            y_pos = dat['/MAPS/y_axis'][:]
            theta = round(float(dat['MAPS/extra_pvs'][1,3].decode('utf-8')))  #get angle int info
            suc = True
    print(f'call_count is {getElmMap.counter}, scan_num is {scan_num}')
    getElmMap.counter += 1
    return getElmMap.counter, suc, elmmap, x_pos, y_pos, theta, f
getElmMap.counter = 0
#----------------------xyl------------------------------------

def kmean_analysis(n_clusters, data, random_state = 52, sigma = None, cval = None,
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
        plt.tight_layout()
        plt.show(block=False)
        
        if (savefig == 1) & (fname is not None):
            fig.savefig(fname, dpi = 300)
    
    return km_label, km_bool, fig


def phaseCorrelate(image1, image2):
    fft_array1 = fft2(image1)
    fft_array2 = fft2(image2)

    shape = image1.shape
    c = abs(ifft2(fft_array1 * fft_array2.conjugate() / (abs(fft_array1) * abs(fft_array2))))
    t0, t1 = np.unravel_index(np.argmax(c), image1.shape)
    if t0 > shape[0] // 2:
        t0 -= shape[0]
    if t1 > shape[1] // 2:
        t1 -= shape[1]
    return t0, t1

# second k-mean analysis - xyl
def kmean_analysis_2(file, elmmap,n_clusters,Gaussian_blur,log_his,sel_cluster,savefolder,theta):  #file can be raw image, log image
    #fig = plt.figure()
    fig, (ax1,ax2,ax3) = plt.subplots(1,3)
    #elmmap = io.imread(elmmap)
    elmmap[np.isnan(elmmap)] = 1e-5
    elmmap[np.isinf(elmmap)] = 1e-5
    h,w = elmmap.shape
    histogram, bin_edges = np.histogram(elmmap, bins=256)
    if log_his==0:    
        img_1d = elmmap.reshape(h*w,1)     
        kmeans = KMeans(n_clusters=n_clusters, init='k-means++',n_init='auto',max_iter=500,tol=1e-4, random_state=0, algorithm='lloyd').fit(img_1d)
        center = kmeans.cluster_centers_
        ax2.plot(bin_edges[0:-1], histogram,linewidth=3)
    else:
        #log histogram        
        log_his = np.log(histogram)
        log_his[np.isnan(log_his)] = 1e-5
        log_his[np.isinf(log_his)] = 1e-5
        bin_log_his = np.vstack((bin_edges[0:-1], log_his)).T
        kmeans = KMeans(n_clusters=n_clusters, init='k-means++',n_init='auto',max_iter=500,tol=1e-4, random_state=0, algorithm='lloyd').fit(bin_log_his)
        center = kmeans.cluster_centers_
        ax2.plot(bin_edges[0:-1], log_his,linewidth=3)
    his_fit = kmeans.fit_predict(img_1d)
    his_fit_2d = his_fit.reshape(h,w)
    if sel_cluster == 0:
        his_fit_2d = 1-his_fit_2d
    else:
        his_fit_2d = his_fit_2d
    if Gaussian_blur != 0:
        his_fit_2d = cv.blur(his_fit_2d,(Gaussian_blur,Gaussian_blur))
    else:
        his_fit_2d = his_fit_2d
    img_cen = ndimage.center_of_mass(his_fit_2d) #(y,x)
    ax1.imshow(elmmap,aspect='auto')
    ax1.plot(img_cen[1], img_cen[0], marker='o',markersize=10)  
    ax2.axvline(center[0,0],color='black',linewidth=2)
    ax2.axvline(center[1,0],color='black',linewidth=2)
    ax3.imshow(his_fit_2d,aspect='auto')
    ax3.plot(img_cen[1], img_cen[0], marker='o',markersize=10)
    fname = os.path.splitext(os.path.basename(file))[0]
    ax1.set_title(file,fontdict={'fontsize':8})
    ax2.set_title(f'log_hist={log_his}_hist+K-meanCenter',fontdict={'fontsize':8})
    ax3.set_title('K-mean_seg',fontdict={'fontsize':8})
    plt.tight_layout()
    #plt.show(block=False)
    fig.savefig(f'{savefolder}/{fname}_ang{theta}.png')
    print('coarse {fname} image saved')  #fname is coarse scan name
    #plt.close(fig)
    return center, his_fit_2d, img_cen

def getROIcoordinate_data_2(fname, elm, n_clusters,Gaussian_blur,log_his,rev_pixval):
    suc, elmmap, x_pos, y_pos = getElmMap(fname, elm)
    if suc == True:
        cen, his_fit_2d, img_cen = kmean_analysis_2(elmmap,n_clusters,Gaussian_blur,log_his,rev_pixval)
        x_cen = int(img_cen[1])
        y_cen = int(img_cen[0]) 
        x_len_range = np.arange(elmmap.shape[1])
        y_len_range = np.arange(elmmap.shape[0])
        #zip_pix_pos = {x_len_range[i]: x_pos[i] for i in range(len(x_len_range))}
        zip_pix_pos_x = dict(zip(x_len_range, x_pos))  
        zip_pix_pos_y = dict(zip(y_len_range, y_pos))  
        new_x = zip_pix_pos_x.get(x_cen)
        new_y= zip_pix_pos_y.get(y_cen)
        width = x_pos[-1] - x_pos[0]
        height = y_pos[-1] - y_pos[0]
        return new_x, new_y, width, height
    else:
        return None
        
    
