from skimage.registration import phase_cross_correlation
import numpy as np
import os
import matplotlib.pyplot as plt
from plotting import xrfOverlay
from pystackreg import StackReg
from scipy.fftpack import fft2, ifft2
from scipy import optimize

def zero_pad(a, shift, axis = None):
    a = np.asanyarray(a)
    if shift == 0: return a
    if axis is None:
        n = a.size
        reshape = True
    else:
        n = a.shape[axis]
        reshape = False
    if np.abs(shift) > n:
        res = np.zeros_like(a)
    elif shift < 0:
        shift += n
        zeros = np.zeros_like(a.take(np.arange(n-shift), axis))
        res = np.concatenate((a.take(np.arange(n-shift,n), axis), zeros), axis)
    else:
        zeros = np.zeros_like(a.take(np.arange(n-shift,n), axis))
        res = np.concatenate((zeros, a.take(np.arange(n-shift), axis)), axis)
    if reshape:
        return res.reshape(a.shape)
    else:
        return res


def elmArrAlignment(elmdata_org, refidx, elms, angles, elm_idx = 0, savefig = False, figsize = (5,5), padzero=True, upsample_factor=1, plotOption = True):
    xshift = np.zeros(angles.shape)
    yshift = np.zeros(angles.shape)
    elmdata = elmdata_org.copy()

    for i in range(elmdata.shape[1]):
        if (i != refidx):
            if i < refidx:
                fidx = refidx-i
                movidx = refidx-i-1
            elif i > refidx:
                fidx = i - 1
                movidx = i

            ref = elmdata[elm_idx,fidx,:,:]
            mov = elmdata[elm_idx,movidx,:,:].copy()
            vmaxref = np.nanmean(ref)+3*np.nanstd(ref)
            vmaxmov = np.nanmean(mov)+3*np.nanstd(mov)

            s = phase_cross_correlation(ref, mov) if upsample_factor == 0 else phase_cross_correlation(ref, mov, upsample_factor = upsample_factor)
            xshift[movidx] = int(s[0][1])
            yshift[movidx] = int(s[0][0])

            for j in range(elmdata.shape[0]):
                mov_ = elmdata[j,movidx,:,:].copy()
                if padzero:
                    reg = zero_pad(zero_pad(mov_.copy(), int(s[0][1]), axis = 1), int(s[0][0]), axis = 0)
                else:
                    reg = np.roll(mov_.copy(), (int(s[0][1]), int(s[0][0])), axis=(1,0))
                    
#                 temp = np.zeros(reg.shape, order = 'c')
#                 sx = abs(int(s[0][1]))
#                 sy = abs(int(s[0][0]))
#                 if (sx > 0) & (sy > 0):
#                     temp[sx:-sx,sy:-sy] = reg[sx:-sx,sy:-sy]
#                 else:
#                     temp = reg
                elmdata[j,movidx,:,:] = reg
            
            if padzero:
                reg = zero_pad(zero_pad(mov.copy(), int(s[0][1]), axis = 1), int(s[0][0]), axis = 0)
            else:
                reg = np.roll(mov.copy(),(int(s[0][1]), int(s[0][0])), axis=(1,0))
                
            if plotOption:
                fig, axes = plt.subplots(2,2, figsize=figsize)
                axes[0,0].imshow(ref,cmap = 'coolwarm_r', vmax = vmaxref)
                axes[0,1].imshow(mov, cmap = 'coolwarm', vmax = vmaxmov)
                axes[0,0].set_xticks([])
                axes[0,0].set_yticks([])
                axes[0,1].set_xticks([])
                axes[0,1].set_yticks([])
                axes[0,0].set_title('Theta = %.2f'%(angles[fidx]))
                axes[0,1].set_title('Theta = %.2f'%(angles[movidx]))
                f = xrfOverlay(ref, mov, vmax_m1 = vmaxref, vmax_m2 = vmaxmov, labels=[str(angles[fidx]),str(angles[movidx])], 
                        cmap=['coolwarm_r','coolwarm'], ax = axes[1,0])
                f = xrfOverlay(ref, reg, vmax_m1 = vmaxref, vmax_m2 = vmaxmov, labels=[str(angles[fidx]),str(angles[movidx])], 
                            cmap=['coolwarm_r','coolwarm'], ax = axes[1,1])
                axes[1,0].set_title('Raw Data')
                axes[1,1].set_title('Aligned Data')
                plt.tight_layout()

                if savefig:
                    img_path= os.path.join(os.getcwd(),'alignment_%s'%(elms[elm_idx]))
                    if not os.path.exists(img_path):
                        os.makedirs(img_path)
                    fpath = os.path.join(img_path, 'theta = %.2f x = %d y = %d.png'%(angles[movidx], xshift[movidx], yshift[movidx]))
                    fig.savefig(fpath, dpi = 300, transparent=True)              
    return elmdata


def elmArrAlignment_pystackreg(elmdata_org, refidx, elms, angles, regType = StackReg.TRANSLATION, elm_idx = 0, savefig = False, figsize = (5,5)):
    xshift = np.zeros(angles.shape)
    yshift = np.zeros(angles.shape)
    elmdata = elmdata_org.copy()
    sr = StackReg(regType)

    for i in range(elmdata.shape[1]):
        if (i != refidx):
            if i < refidx:
                fidx = refidx-i
                movidx = refidx-i-1
            elif i > refidx:
                fidx = i - 1
                movidx = i

            ref = elmdata[elm_idx,fidx,:,:]
            mov = elmdata[elm_idx,movidx,:,:].copy()
            vmaxref = np.nanmean(ref)+3*np.nanstd(ref)
            vmaxmov = np.nanmean(mov)+3*np.nanstd(mov)

            tform = sr.register(ref, mov)
            for j in range(elmdata.shape[0]):
                mov_ = elmdata[j,movidx,:,:].copy()
                reg = sr.transform(mov_)
                elmdata[j,movidx,:,:] = reg
            
            reg = sr.transform(mov)
                
            fig, axes = plt.subplots(2,2, figsize=figsize)
            axes[0,0].imshow(ref,cmap = 'coolwarm_r', vmax = vmaxref)
            axes[0,1].imshow(mov, cmap = 'coolwarm', vmax = vmaxmov)
            axes[0,0].set_xticks([])
            axes[0,0].set_yticks([])
            axes[0,1].set_xticks([])
            axes[0,1].set_yticks([])
            axes[0,0].set_title('Theta = %.2f'%(angles[fidx]))
            axes[0,1].set_title('Theta = %.2f'%(angles[movidx]))
            f = xrfOverlay(ref, mov, vmax_m1 = vmaxref, vmax_m2 = vmaxmov, labels=[str(angles[fidx]),str(angles[movidx])], 
                    cmap=['coolwarm_r','coolwarm'], ax = axes[1,0])
            f = xrfOverlay(ref, reg, vmax_m1 = vmaxref, vmax_m2 = vmaxmov, labels=[str(angles[fidx]),str(angles[movidx])], 
                        cmap=['coolwarm_r','coolwarm'], ax = axes[1,1])
            axes[1,0].set_title('Raw Data')
            axes[1,1].set_title('Aligned Data')
            plt.tight_layout()
            
            if savefig:
                img_path= os.path.join(os.getcwd(),'alignment_%s'%(elms[elm_idx]))
                if not os.path.exists(img_path):
                    os.makedirs(img_path)
                fpath = os.path.join(img_path, 'theta = %.2f x = %d y = %d.png'%(angles[movidx], xshift[movidx], yshift[movidx]))
                fig.savefig(fpath, dpi = 300, transparent=True)              
    return elmdata


def checkAlignment(elmdata, refidx, elms, angles, elm_idx = 0, figsize = (5,5), savefig=False):
    for i in range(elmdata.shape[1]):
        if (i != refidx):
            if i < refidx:
                fidx = refidx-i
                movidx = refidx-i-1
            elif i > refidx:
                fidx = i - 1
                movidx = i

            ref = elmdata[elm_idx,fidx,:,:]
            mov = elmdata[elm_idx,movidx,:,:].copy()
            vmaxref = np.nanmean(ref)+3*np.nanstd(ref)
            vmaxmov = np.nanmean(mov)+3*np.nanstd(mov)

            fig, axes = plt.subplots(1,3, figsize=figsize)
            axes[0].imshow(ref,cmap = 'coolwarm_r', vmax = vmaxref)
            axes[1].imshow(mov, cmap = 'coolwarm', vmax = vmaxmov)
            axes[0].set_xticks([])
            axes[0].set_yticks([])
            axes[1].set_xticks([])
            axes[1].set_yticks([])
            axes[0].set_title('Theta = %.2f'%(angles[fidx]))
            axes[1].set_title('Theta = %.2f'%(angles[movidx]))
            f = xrfOverlay(ref, mov, vmax_m1 = vmaxref, vmax_m2 = vmaxmov, labels=[str(angles[fidx]),str(angles[movidx])], 
                    cmap=['coolwarm_r','coolwarm'], ax = axes[2])
            axes[2].set_title('Overlay Data')
            plt.tight_layout()
            
            if savefig:
                img_path= os.path.join(os.getcwd(),'overlay_%s'%(elms[elm_idx]))
                if not os.path.exists(img_path):
                    os.makedirs(img_path)
                fpath = os.path.join(img_path, 'theta = %.2f x = %d y = %d.png'%(angles[movidx]))
                fig.savefig(fpath, dpi = 300, transparent=True)    
                
def xCor(xcorElementIndex, projections, data, xshift, yshift, refscan = None):
    
    for i in np.arange(projections - 1):
        # onlyfilename=self.fileNames[i+1].rfind("/")
        if refscan is None:
            img1 = data[xcorElementIndex, i, :, :]
        else:
            img1 = data[xcorElementIndex, refscan, :, :]
            
        img2 = data[xcorElementIndex, i + 1, :, :]

        t0, t1 = crossCorrelate(img1, img2)
        data[:, i + 1, :, :] = np.roll(data[:, i + 1, :, :], t0, axis=1)
        data[:, i + 1, :, :] = np.roll(data[:, i + 1, :, :], t1, axis=2)
        xshift[i + 1] += t1
        yshift[i + 1] += t0
        
def crossCorrelate(image1, image2):
    '''

        :param image1: 2d array
        :param image2: 2d array
        :return:
    '''
    fft_array1 = fft2(image1)
    fft_array2 = fft2(image2)

    shape = image1.shape
    c = abs(ifft2(fft_array1 * fft_array2.conjugate()))
    t0, t1 = np.unravel_index(np.argmax(c), image1.shape)
    if t0 > shape[0] // 2:
        t0 -= shape[0]
    if t1 > shape[1] // 2:
        t1 -= shape[1]

    return t0, t1