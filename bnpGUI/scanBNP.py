'''
Control PV based on type of measurements

'''
#!/home/beams/USERBNP/.conda/envs/py36/bin/python

import os, time, sys
import numpy as np
from imgProcessing import getROIcoordinate_data, getElmMap, kmean_analysis_2
import matplotlib.pyplot as plt
from skimage import io
from pvComm import pvCommsubclass
from epics import PV, caget, caput, cainfo
def parmLabelToPVdict():
    d = {'width':'x_width', 'height':'y_width', 'w_step':'x_step',
         'h_step':'y_step', 'dwell':'dwell', 'x_scan':'x_center_Rqs', 
         'y_scan':'y_center_Rqs', 'z_scan': 'z_value_Rqs', 'target_theta':'sm_rot_Rqs','Dwell (ms)':'dwell_step'}
    return d

def xrfSetup(pvComm, scandic):
    d = parmLabelToPVdict()  #list d to dictionary d
    #d = {'width': 'x_width',
    #     'height': 'y_width',
    #     'w_step': 'x_step',
    #     'h_step': 'y_step',
    #     'dwell': 'dwell',
    #     'x_scan': 'x_center_Rqs',
    #     'y_scan': 'y_center_Rqs',
    #     'z_scan': 'z_value_Rqs',
    #     'target_theta': 'sm_rot_Rqs'}      key is the name in GUI, value is the name in pvs (pvobjects)
    parms = ['width', 'height', 'w_step', 'h_step', 
             'dwell', 'y_scan', 'x_scan']  #names used in GUI
    parm_label = [d[s] for s in parms]
    #parm_label = ['x_width',
    #              'y_width',
    #              'x_step',
    #              'y_step',
    #              'dwell',
    #              'y_center_Rqs',
    #              'x_center_Rqs'] name in pvs
    parm_value = [float(scandic[s]) for s in parms]  #get value of parameters in parms
    pvComm.writeScanInit('XRF', scandic['smpName'], str(scandic)) #log information
    pvComm.blockBeamBDA(scandic['bda']) #move bda to block position
    # pvComm.changeXYcombinedMode()
    pvComm.changeXtoCombinedMode()  #change x to combine mode
    pvComm.assignPosValToPVs(parm_label, parm_value)  #give value to pvs
    return getMotorList(scandic)
#-----------------------------xyl: xanes setup---------------------------
def xrfSetup_fromxanes(pvComm, scandic):
    pvComm.assignEng('mono_eng', scandic['Energe (keV)'])   # use the energy setting in XANES part
    v_m_e = ['read_1','drive_1','mono_mode','collect_mode']             
    pvComm.assignPosVaToPvs(v_m_e,['9idbXMAP:scan1.R1PV','9idbXMAP:scan1.P1PV',1,1])
    pvComm.assignEng('mono_eng', scandic['Energe (keV)'])   # use the energy setting in XANES part
    
    d = parmLabelToPVdict()  #list d to dictionary d
    #d = {'width': 'x_width',
    #     'height': 'y_width',
    #     'w_step': 'x_step',
    #     'h_step': 'y_step',
    #     'dwell': 'dwell',
    #     'x_scan': 'x_center_Rqs',
    #     'y_scan': 'y_center_Rqs',
    #     'z_scan': 'z_value_Rqs',
    #     'target_theta': 'sm_rot_Rqs'}      key is the name in GUI, value is the name in pvs (pvobjects)
    parms = ['width', 'height', 'w_step', 'h_step', 
             'dwell', 'y_scan', 'x_scan','xanes_eng_cen']  #names used in GUI
    parm_label = [d[s] for s in parms]
    #parm_label = ['x_width',
    #              'y_width',
    #              'x_step',
    #              'y_step',
    #              'dwell',
    #              'y_center_Rqs',
    #              'x_center_Rqs'] name in pvs
    parm_value = [float(scandic[s]) for s in parms]  #get value of parameters in parms
    pvComm.writeScanInit('XRF', scandic['smpName'], str(scandic)) #log information
    pvComm.blockBeamBDA(scandic['bda']) #move bda to block position
    # pvComm.changeXYcombinedMode()
    pvComm.changeXtoCombinedMode()  #change x to combine mode
    pvComm.assignPosValToPVs(parm_label, parm_value)  #give value to pvs
    return getMotorList(scandic)
def xanes_ps_n(pvComm,scandic):   #setup for xanes from xrf:   
    v_m_e = ['read_1','drive_1','mono_mode','collect_mode']             
    pvComm.assignPosVaToPvs(v_m_e,['','9idb:mono_pid1.FBON',0,0])
    pvComm.assignEng('mono_eng', scandic['Energe (keV)'])
    #'mono_mode': '9idb:mono_pid1.FBON', 'read_1':'9idbXMAP:scan1.R1PV', 
    #'drive_1':'9idbXMAP:scan1.P1PV', 'mono_eng':'2ida2:BraggEAO.VAL',
    #'dwell_step': '9idbXMAP:userTran1.P', 'xanes_eng_cen':'9idbXMAP:scan1.P1CP', 'collect_mode':'9idbXMAP:CollectMode',
    
    d = parmLabelToPVdict()  #list d to dictionary 
    #d = {'width': 'x_width',
    #     'height': 'y_width',
    #     'w_step': 'x_step',
    #     'h_step': 'y_step',
    #     'dwell': 'dwell',
    #     'x_scan': 'x_center_Rqs',
    #     'y_scan': 'y_center_Rqs',
    #     'z_scan': 'z_value_Rqs',
    #     'Dwell (ms)':'dwell_step'
    #     'target_theta': 'sm_rot_Rqs'}      key is the name in GUI, value is the name in pvs (pvobjects)
    # width-->'9idbBNP:scan1.P1WD', w_step--> ''9idbBNP:scan1.P1SI'
    parms = ['y_scan', 'x_scan']  #names used in GUI
    parm_label = [d[s] for s in parms] + ['xanes_eng_cen','width', 'w_step','dwell_step']
    #parm_label = ['x_width',
    #              'x_step',
    #              'dwell',
    #              'y_center_Rqs',
    #              'x_center_Rqs',] name in pvs
    name_gui_eng = ['Energy (keV)','Energy width (keV)', 'Energy step (keV)', 'Dwell (s)'] 
    parm_value = [float(scandic[s]) for s in parms] + [float(scandic[n]) for n in name_gui_eng]  #get value of parameters in parms and 0,0 for mono, collect mode
    pvComm.writeScanInit('XANES (fixed region)', scandic['smpName'], str(scandic)) #log information
    pvComm.blockBeamBDA(scandic['bda']) #move bda to block position
    pvComm.changeXtoCombinedMode()
    pvComm.assignPosValToPVs(parm_label, parm_value)  #give value to pvs
    return getMotorList(scandic)

def xanes_ps_c(pvComm,scandic):   #setup for xanes from xrf:    
    d = parmLabelToPVdict()  #list d to dictionary 
    #d = {'width': 'x_width',
    #     'height': 'y_width',
    #     'w_step': 'x_step',
    #     'h_step': 'y_step',
    #     'dwell': 'dwell',
    #     'x_scan': 'x_center_Rqs',
    #     'y_scan': 'y_center_Rqs',
    #     'z_scan': 'z_value_Rqs',
    #     'Dwell (ms)':'dwell_step'
    #     'target_theta': 'sm_rot_Rqs'}      key is the name in GUI, value is the name in pvs (pvobjects)
    # width-->'9idbBNP:scan1.P1WD', w_step--> ''9idbBNP:scan1.P1SI'
    parms = ['y_scan', 'x_scan']  #names used in GUI
    parm_label = [d[s] for s in parms] + ['xanes_eng_cen','width', 'w_step','dwell_step'] #'dwell_step': '9idbXMAP:userTran1.P', 'xanes_eng_cen':'9idbXMAP:scan1.P1CP', 'collect_mode':'9idbXMAP:CollectMode',
    #parm_label = ['x_width',
    #              'x_step',
    #              'dwell',
    #              'y_center_Rqs',
    #              'x_center_Rqs',] name in pvs
    name_gui_eng = ['Energy (keV)','Energy width (keV)', 'Energy step (keV)', 'Dwell (s)'] 
    parm_value = [float(scandic[s]) for s in parms] + [float(scandic[n]) for n in name_gui_eng]  #get value of parameters in parms and 0,0 for mono, collect mode
    pvComm.writeScanInit('XANES (fixed region)', scandic['smpName'], str(scandic)) #log information
    pvComm.blockBeamBDA(scandic['bda']) #move bda to block position
    pvComm.changeXtoCombinedMode()
    pvComm.assignPosValToPVs(parm_label, parm_value)  #give value to pvs
    return getMotorList(scandic)
#-----------------------------------------------------------------------------------------------    
def getMotorList(scandic):
    p = ['target_theta', 'z_scan', 'y_scan', 'x_scan']
    motorlabel = ['sm_rot', 'z_value', 'y_center', 'x_center']
    mtolerance = [0.1, 0.5, 0.1, 0.1]
    mlist = []
    for p_, ml_, mt_ in zip(p, motorlabel, mtolerance):
        mlist.append((ml_, float(scandic[p_]), mt_))
    return mlist


def scanStart(scandic, pvComm, bda):   #revised added xanes where not changing to pixel mode
    if scandic['scanType'] == 'XANES (fixed region)':
        pvComm.openBeamBDA(bda)
    else:
        pvComm.changeXtoPiezolMode()
        pvComm.openBeamBDA(bda)

def scanFinish(scandic, pvComm, bda):
    if scandic['scanType'] == 'XANES (fixed region)':
        time.sleep(3)
        pvComm.blockBeamBDA(bda)
        pvComm.changeXtoCombine()
        time.sleep(0.5)
    else:
        time.sleep(5)
        pvComm.blockBeamBDA(bda)
        pvComm.changeXtoCombinedMode()
        pvComm.centerPiezoY(waittime=3)   #at 2idd, stage issue, work around
        pvComm.centerPiezoY(waittime=3)   #at 2idd, stage issue, work around
        time.sleep(0.5)
    
def fileReady(coarse_sc, fdir, tlim = 30):
    fpath = os.path.join(fdir, 'img.dat/%s.h5'%(coarse_sc))
    if os.path.exists(fpath):
        fmtime = os.path.getmtime(fpath)
        tdiff = time.time() - fmtime
        if tdiff > tlim:
            return 1
        else:
            sys.stdout.write('Waiting for coarse scan file %s.h5 to be ready,'\
                     ' file modified time: %d, time difference: %d \n'\
                     %(coarse_sc, fmtime, tdiff))
            return 0
    else:
        sys.stdout.write('File %s not exisit\n'%fpath)
        return 0
    
def imgProgFolderCheck(fdir):
    img_path= os.path.join(fdir,'imgProg')
    if not os.path.exists(img_path):
        os.makedirs(img_path)
    return img_path
'''
def getCoordinate(pvComm, coarse_sc, scandic, n_std = 2):
    fready = fileReady(coarse_sc, pvComm.userdir)
    coarse_h5path = os.path.join(pvComm.userdir, 'img.dat/%s.h5'%(coarse_sc))
    if fready:
        imgfolder = imgProgFolderCheck(pvComm.userdir)   
        imgpath = os.path.join(imgfolder, 'bbox_%s.png'%(coarse_sc))
        print(imgpath)
        elmmap = getElmMap(coarse_h5path, scandic['elm'])
        
        mask = np.ones(elmmap[0].shape)
        if scandic['use_mask']:
            maskmap = getElmMap(coarse_h5path, scandic['mask_elm'])
            mask = maskmap < (np.mean(maskmap) + n_std * np.std(maskmap.ravel()))
        m = elmmap[0] * mask
        
        x, y, w, h = getROIcoordinate_data(m, elmmap[1], elmmap[2], 
                                           n_cluster = scandic['n_clusters'],
                                           sel_cluster = scandic['sel_cluster'],
                                           figpath = imgpath)
        return np.round(x,2), np.round(y,2)
    else:
        return None
'''   
#------------------------------------------------xyl----------------
def getCoordinate_2(coarse_sc, scandic,flag):
    uf = pvCommsubclass().user_di()
    coarse_h5 = os.path.join(uf,'img.dat/',coarse_sc)
    
    
    #fready = fileReady(coarse_sc, pvComm.userdir)
    #coarse_h5path = os.path.join(pvComm.userdir, 'img.dat/%s.h5'%(coarse_sc))  
    if os.path.exists(coarse_h5):   
        if flag == 0:
            count, suc, elmmap, x_pos, y_pos,theta, f = getElmMap(fname=coarse_h5, elm=scandic['elm'],flag=0)
            #fready = fileReady(coarse_sc, pvComm.userdir)
            #coarse_h5path = os.path.join(pvComm.userdir, 'img.dat/%s.h5'%(coarse_sc))       
            #if fready:
            #imgfolder = imgProgFolderCheck(pvComm.userdir)   
            #imgpath = os.path.join(imgfolder, 'bbox_%s.png'%(coarse_sc))
           # print(imgpath)
            #count, suc, elmmap, x_pos, y_pos,f = getElmMap(fname=coarse_h5path, elm=scandic['elm'],flag=flag)
            #imgfolder = imgProgFolderCheck(pvComm.userdir)
            #file_name = os.path.splitext(os.path.basename(coarse_sc))[0]
            #io.imsave(f'{savefolder}/{file_name}.tif',elmmap)   #try to save image in save folder
        #else:
          #  return None
        else:
            test_sc_range = np.arange(55,160,2)   #for test
            test_folder = '/mnt/micdata1/bnp/2022-3/Merk/img.dat.5det.ele/'
            count, suc, elmmap, x_pos, y_pos,theta, f = getElmMap(elm=scandic['elm'], test_folder=test_folder, test_sc_range=test_sc_range, flag=1)

        if suc == True:
            try:
                cen, his_fit_2d, img_cen = kmean_analysis_2(file=coarse_h5,savefolder=f'{pvCommsubclass().user_folder()}', elmmap=elmmap,n_clusters=scandic['n_clusters'],Gaussian_blur=scandic['Gaussian_blur'],
                                                        log_his=scandic['log_his'],sel_cluster=scandic['sel_cluster'],theta=theta)
                x_cen = round(img_cen[1])
                y_cen = round(img_cen[0])                              
                x_len_range = np.arange(elmmap.shape[1])
                y_len_range = np.arange(elmmap.shape[0])
            #zip_pix_pos = {x_len_range[i]: x_pos[i] for i in range(len(x_len_range))}
                zip_pix_pos_x = dict(zip(x_len_range, x_pos))  
                zip_pix_pos_y = dict(zip(y_len_range, y_pos))  
                new_x = zip_pix_pos_x.get(x_cen)
                new_y= zip_pix_pos_y.get(y_cen)
            #width = x_pos[-1] - x_pos[0]
            #height = y_pos[-1] - y_pos[0]
                return np.round(new_x,2), np.round(new_y,2),theta  #theta only for test images
            except ValueError:
                return None
        else:
            return None
    else:
        return None

    
        
    
            




