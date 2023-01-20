'''
This is the scan class for data collection at BNP. 

TO DO: add temp PV and output to logbook; remove unimportant lines in logbook
'''

import time, os, sys
from tqdm import tqdm
from epics import caput, caget
from time import gmtime, strftime
import numpy as np
from imgProcessing import *
import pandas as pd

class bnpScan():
    
    def __init__(self, userdir, logfile):
        self.userdir = userdir
        self.logfilepath = os.path.join(userdir, logfile)
        self.logfid = open(self.logfilepath, 'a')
        self.pvs = self.definePVs()
        self.scandic = None
    
    def logger(self, msg):
        sys.stdout.write(msg)
        sys.stdout.flush()
        if self.logfid.closed:
            self.logfid = open(self.logfilepath, 'a')
        self.logfid.write(msg)
        self.logfid.flush()
        
    def getCurrentTime(self):
        ts = pd.Timestamp.now()
        ts_str = ts.strftime('%Y-%m-%d %X')
        return ts_str
    
    def checkPtsArea(self, pts_area, label):
        if len(pts_area) != 5:
            raise ValueError('%s input is invalid. pts_area'\
                             '(x-width, y-width, x-step, y-step, dwell(ms))\n'%(label))
        
    def checkOrgPos_0theta(self, orgPos_xyz_0theta, label):
        if len(orgPos_xyz_0theta) != 3:
            raise ValueError('%s input is invalid. '\
                             '%s (x-center, y-center, z-center) at theta 0\n'%(label, label))
    
    def setupAngleSweepScans(self, sampleName, scans, pts_area, orgPos_xyz_0theta, BDAin):
        # check if inputs are valid
        for s_ in scans:
            if (isinstance(s_, float)) | (isinstance(s_, int)):
                pass
            else:
                raise ValueError('Check input parameters for scans. It should be a list of angles\n')
        
        self.checkPtsArea(pts_area, 'pts_area')
        self.checkOrgPos_0theta(orgPos_xyz_0theta, 'orgPos_xyz_0theta')
        scandic = {'scanMode': 'angleSweep', 'BDAin':BDAin, 'elm':None, 'smpInfo':sampleName,
                   'angleSweep':{'pts_area':pts_area, 'orgPos_xyz_0theta':orgPos_xyz_0theta,
                                 'scans': scans, 'pre_parm':['pts_area', 'orgPos_xyz_0theta'],
                                 'parm_label':['x_width', 'y_width', 'x_step', 'y_step', 'dwell',
                                               'x_center_Rqs', 'y_center_Rqs', 'z_value_Rqs',
                                               'tomo_rot_Rqs']}}
        self.scandic = scandic
    
    def setupBatchXRFScans(self, sampleName, scans, BDAin, smp_angle = None):
        err_msg = 'Input parameters for scans is invalid. Scans input is a list of list.\n'\
                  'Scans: [[x-width, y-width, x-step, y-step, dwell (ms), x-center,'\
                  ' y-center, z-position]]\n'
        for s_ in scans:
            if (isinstance(s_, list)):
                if (len(s_) == 8): pass
                else: raise ValueError(err_msg)
            else:
                raise ValueError(err_msg)
                
        scandic = {'scanMode':'batchXRF_fixAngle', 'BDAin':BDAin, 'elm':None, 'smpInfo':sampleName, 
                  'batchXRF_fixAngle':{'scans': scans, 'smp_angle':smp_angle,
                                       'pre_parm':None,
                                       'parm_label':['x_width', 'y_width', 'x_step', 'y_step',
                                                     'dwell', 'x_center_Rqs', 'y_center_Rqs',
                                                     'z_value_Rqs']}}
        self.scandic = scandic
        
    def setupCoarseFineScans(self, sampleName, scans, BDAin, pts_area_coarse, orgPos_xyz_0theta, 
                             elm, pts_area_fine, n_cluster = 2, sel_cluster = 1, fine_bbox = True):
        self.checkPtsArea(pts_area_coarse, 'pts_area_coarse')
        self.checkPtsArea(pts_area_fine, 'pts_area_fine')
        self.checkOrgPos_0theta(orgPos_xyz_0theta, 'orgPos_xyz_0theta')
        for s_ in scans:
            if (isinstance(s_, float)) | (isinstance(s_, int)):
                pass
            else:
                raise ValueError('Check input parameters for scans. It should be a list of angles\n')
                
        scandic = {'scanMode':'coarse_fine', 'elm':elm, 'smpInfo':sampleName, 'BDAin': BDAin,
                   'coarse_fine':{'pts_area':pts_area_coarse, 'orgPos_xyz_0theta':orgPos_xyz_0theta,
                       'scans':scans, 'pre_parm':['pts_area', 'orgPos_xyz_0theta'],
                       'parm_label':['x_width', 'y_width', 'x_step', 'y_step', 'dwell', 
                                    'x_center_Rqs','y_center_Rqs', 'z_value_Rqs', 'tomo_rot_Rqs'],
                       'find_bbox':True, 'fine_pts_area':pts_area_fine,
                       'n_cluster':n_cluster, 'sel_cluster':sel_cluster}}
        self.scandic = scandic
    
    def definePVs(self):
        self.logger('\n\n%s: Associate motors with PVs\n'%(self.getCurrentTime()))
        pvs = {'x_center_Rqs':'9idbTAU:SM:PX:RqsPos', 'x_center_Act':'9idbTAU:SM:PX:ActPos',
               'y_center_Rqs':'9idbTAU:SY:PY:RqsPos', 'y_center_Act':'9idbTAU:SY:PY:ActPos',
               'z_value_Rqs':'9idbTAU:SM:SZ:RqsPos', 'z_value_Act':'9idbTAU:SM:SZ:ActPos',
               'tomo_rot_Rqs':'9idbTAU:SM:CT:RqsPos', 'tomo_rot_Act':'9idbTAU:SM:CT:ActPos',
               'sm_rot_Rqs':'9idbTAU:SM:ST:RqsPos', 'sm_rot_Act':'9idbTAU:SM:ST:ActPos',
               'x_width':'9idbBNP:scan1.P1WD', 'y_width':'9idbBNP:scan2.P1WD',
               'x_step':'9idbBNP:scan1.P1SI', 'y_step':'9idbBNP:scan2.P1SI',
               'dwell':'9idbBNP:scanTran3.C', 'BDA_pos':'9idbTAU:UA:UX:RqsPos',
               
               'x_motorMode':'9idbTAU:SM:Ps:xMotionChoice.VAL',
               'y_motorMode':'9idbTAU:SY:Ps:yMotionChoice.VAL',
               'x_setcenter':'9idbBNP:aoRecord11.PROC', 'y_setcenter':'9idbBNP:aoRecord12.PROC',
               'piezo_xCenter':'9idbTAU:SM:Ps:xCenter.PROC',
               'piezo_yCenter':'9idbTAU:SY:Ps:yCenter.PROC',
               'tot_lines':'9idbBNP:scan2.NPTS', 'cur_lines':'9idbBNP:scan2.CPT',
               'temp':'9idbCRYO:CryoCon1:In_3:Temp.VAL',
               
               'CryoCon1:In_1':'9idbCRYO:CryoCon1:In_1:Temp.VAL',
               'CryoCon1:In_3':'9idbCRYO:CryoCon1:In_3:Temp.VAL',
               'CryoCon1:In_2':'9idbCRYO:CryoCon1:In_2:Temp.VAL',
               'CryoCon3:In_2':'9idbCRYO:CryoCon3:In_2:Temp.VAL',
               'CryoCon3:Loop_2':'9idbCRYO:CryoCon3:Loop_2:SetControl.VAL',
               
               
               'run':'9idbBNP:scan2.EXSC',
               'fname_saveData':'9idbBNP:saveData_fileName',
               'next_scnum':'9idbBNP:saveData_scanNumber'
               }
        return pvs
    
    def changeTomoRotate(self, theta):
        curr_angle = caget(self.pvs['tomo_rot_Act'])
        t = self.getCurrentTime()
        self.logger('%s; Changing tomo rotation angle from to %.2f to %.2f\n'%(t, curr_angle, theta))
        caput(self.pvs['tomo_rot_Rqs'], theta)
        self.motorReady(['tomo_rot'],[0.1])
    
    def changeSMRotate(self, theta):
        curr_angle = caget(self.pvs['sm_rot_Act'])
        t = self.getCurrentTime()
        self.logger('%s; Changing sample rotation angle from to %.2f to %.2f\n'%(t, curr_angle, theta))
        caput(self.pvs['sm_rot_Rqs'], theta)
        self.motorReady(['sm_rot'],[0.1])
    
    def blockBeamBDA(self):
        bda_pos = self.scandic['BDAin'] - 500
        t = self.getCurrentTime()
        self.logger('%s: Move BDA to block position at: %.3f\n'%(t, bda_pos))
        caput(self.pvs['BDA_pos'], bda_pos)
        
    def openBeamBDA(self):
        bda_pos = self.scandic['BDAin']
        self.logger('%s: Move BDA to open position at: %.3f\n'%(self.getCurrentTime(), bda_pos))
        caput(self.pvs['BDA_pos'], bda_pos)
    
    def changeXYcombinedMode(self):        
        self.logger('%s; Changing XY scan mode to combined motion\n'%(self.getCurrentTime()))
        caput(self.pvs['x_motorMode'], 0)
        time.sleep(2.)
        caput(self.pvs['y_motorMode'], 0)
        time.sleep(2.)
        
    def changeXtoPiezolMode(self):
        self.logger('%s: Changing X scan mode to Piezo only\n'%(self.getCurrentTime()))
        caput(self.pvs['x_motorMode'], 2)
        time.sleep(1.)

    def setXYcenter(self):
        self.logger('%s: Update the current position as the center of'\
                    'the scan.\n'%(self.getCurrentTime()))
        caput(self.pvs['x_setcenter'], 1)
        caput(self.pvs['y_setcenter'], 1)
        time.sleep(0.1)
    
    def assignPosValToPVs(self, pvstr, pvval):
        for s_, v_ in zip(pvstr, pvval):
            caput(self.pvs[s_], v_)
            time.sleep(1)
            self.logger('%s: Change %s to %.3f\n' % (self.getCurrentTime(), s_, v_))

    def centerPiezoXY(self):
        self.logger('%s: Centering piezoX and piezoY.\n'%(self.getCurrentTime()))
        for i in range(2):
            caput(self.pvs['piezo_xCenter'], 1)
            time.sleep(1.)
            caput(self.pvs['piezo_yCenter'], 1)
            time.sleep(1.)
            
    def motorReady(self, label, mtolerance):
        self.logger('%s: Checking whether motors are ready.\n'%(self.getCurrentTime()))
        ready = 0
        rules = [0] * len(label)
        while not ready:
            self.logger('%s: Motors not in position\n'%(self.getCurrentTime()))
            for i, l_ in enumerate(label):
                caput(self.pvs['%s_Rqs'%(l_)], caget(self.pvs['%s_Rqs'%(l_)]))
                time.sleep(1)
                pos_diff = abs(caget(self.pvs['%s_Rqs'%(l_)]) - caget(self.pvs['%s_Act'%(l_)]))
                rules[i] = (pos_diff < mtolerance[i])
                if rules[i]:
                    self.logger('%s: %s motor is in position with value'\
                                '%.2f um\n'%(self.getCurrentTime(), l_, caget(self.pvs['%s_Rqs'%(l_)])))
            
            if all(rules):
                ready = 1
                    
        self.logger('%s: Motors Ready \n'%(self.getCurrentTime()))
        
    def fileReady(self, next_sc, filepath, time_lim):
        # Wait until file exists
        while not os.path.exists(filepath):
            time.sleep(1)

        time_diff = 0
        while (time_diff < time_lim):
            time.sleep(1)
            file_mod_time = os.stat(filepath).st_mtime
            time_diff = int(time.time() - file_mod_time)
            sys.stdout.write('Waiting for coarse scan file %s to be ready,'\
                             ' file modified time: %d, time difference: %d \n'\
                             %(next_sc, file_mod_time, time_diff))
            
    def imgProgFolderCheck(self):
        img_path= os.path.join(self.userdir,'imgProg')
        if not os.path.exists(img_path):
            os.makedirs(img_path)
        return img_path

    def execScan(self, scname, scidx, n_scns):
        self.setXYcenter()
        self.changeXtoPiezolMode()
        self.centerPiezoXY()
        self.openBeamBDA()
        
        # Execute Scan
        caput(self.pvs['run'], 1)
        self.logger('%s: Scanning \n'%(self.getCurrentTime()))
        nlines = caget(self.pvs['tot_lines'])
        tic = time.perf_counter()
        tic1 = tic
        logpvs = ['CryoCon1:In_3', 'CryoCon1:In_2', 'CryoCon3:In_2', 'CryoCon3:Loop_2', 'CryoCon1:In_1']

        with tqdm(total = nlines, desc = '%s (%d/%d)'%(scname, scidx, n_scns)) as pbar:
            pbar.update(caget(self.pvs['cur_lines']))
            
            while caget(self.pvs['run']):
                toc = time.perf_counter()
                if (toc-tic1) >= 1:
                    msg = self.getCurrentTime() + ': '
                    cline = caget(self.pvs['cur_lines'])
                    msg += '%d, '%(cline)
                    for lpv in logpvs:
                        msg += '%.3f, '%(caget(self.pvs[lpv]))
                    msg += '\n'
                    self.logger(msg)
                    tic1 = toc

                elif (toc-tic) >=10:
    #             time.sleep(10)
                    cline = caget(self.pvs['cur_lines'])
                    pbar.update(cline-pbar.n)
                    sys.stdout.write('Scanning %s (batch %d/%d): line %d/%d is done\n'\
                                     %(scname, scidx+1, n_scns, cline, nlines))
                    tic = toc

        self.logger('%s: Finish scan: %s%s'%(self.getCurrentTime(), scname, '\n'*3))
        self.blockBeamBDA()

    def fineScanInit(self, params_label, params, scname, scidx, n_scns):
        self.blockBeamBDA()
        self.changeXYcombinedMode()
        self.assignPosValToPVs(params_label, params)

        # Check following motor to see if they are in position
        label = ['x_center', 'y_center', 'z_value', 'sm_rot']
        mtolerance = [0.1, 0.1, 0.5, 0.2]
        self.motorReady(label, mtolerance)
        self.execScan(scname, scidx, n_scns)
        
    def batchXRFInit(self, params_label, params, scname, scidx, n_scns):
        self.blockBeamBDA()
        self.changeXYcombinedMode()
        
        if self.scandic['batchXRF_fixAngle']['smp_angle'] is not None:
            self.changeSMRotate(self.scandic['batchXRF_fixAngle']['smp_angle'])
            
        self.assignPosValToPVs(params_label, params)
        label = ['x_center', 'y_center', 'z_value']
        mtolerance = [0.1, 0.1, 0.5]
        self.motorReady(label, mtolerance)
        self.execScan(scname, scidx, n_scns)
    
    def angleSweepScanInit(self, params_label, params, scname, scidx, n_scns):
        self.blockBeamBDA()
        self.changeXYcombinedMode()
        
        self.logger('%s: Putting Tomo angle rotation back to 0 deg\n'%(self.getCurrentTime()))
        self.changeTomoRotate(0) # initialize angle at 0 deg
        time.sleep(1)
        
        self.assignPosValToPVs(params_label[:-1], params[:-1])
        self.motorReady(['x_center', 'y_center', 'z_value'], [0.1, 0.1, 0.5])
        self.changeTomoRotate(params[-1]) # change to desired angle
        self.motorReady(['x_center', 'z_value'], [0.1, 0.5])
        self.execScan(scname, scidx, n_scns)
        
    def coarseFineScanInit(self, params_label, params, scname, scidx, n_scns, scan_setting):
        self.angleSweepScanInit(params_label, params, scname, scidx, n_scns)
        
        # If find_bbox is on, perform image processing on the coarse scan to get coordinates
        if scan_setting['find_bbox']:
            cscan_path = os.path.join(self.userdir, 'img.dat/%s.h5'%(scname))
            time_lim = 10  #sec
            self.fileReady(scname, cscan_path, time_lim)
            img_path = self.imgProgFolderCheck()   
            figpath = os.path.join(img_path, 'bbox_%s.png'%(scname))
            new_x, new_y, new_w, new_h = getROIcoordinate(cscan_path, 
                                                          self.scandic['elm'],
                                                          n_cluster = scan_setting['n_cluster'],
                                                          sel_cluster = scan_setting['sel_cluster'],
                                                          figpath = figpath)
            f_scanparm = scan_setting['fine_pts_area']
            
#             ##TODO: implement checkROIIntensity function
#             proceed = checkROIIntensity(cscan_path, self.scandic['elm'], 
#                                         [new_x, new_y, f_scanparm[0], f_scanparm[1]])
            proceed = 1
            
            if proceed:
                self.logger('%.2f(width) \n %.2f(height)\n'%(new_w, new_h))
                curr_smz = caget(self.pvs['z_value_Rqs'])
                scan_ = params[-1]
                params = []
                f_scanparm = scan_setting['fine_pts_area']
                params = f_scanparm + [scan_] + [new_x, new_y, curr_smz]
                flabels = ['x_width', 'y_width', 'x_step', 'y_step', 'dwell', 'sm_rot_Rqs',
                              'x_center_Rqs','y_center_Rqs', 'z_value_Rqs']

                next_sc = self.nextScanName(caget(self.pvs['fname_saveData']))
                self.logger('%s Initiating fine scan %s %s\n'%('#'*20, next_sc, '#'*20))
                self.logger('Sample temp (K): %.3f\n'%(caget(self.pvs['temp'])))
                self.fineScanInit(flabels, params, next_sc, scidx, n_scns)

            else:
                self.logger('Extracted ROI appears to have intensity below average,'\
                            ' suggesting a no feature region.\n Aborting the batch scan. \n')
                ##TODO: send out an email if this happen
                status = -1  # Flag to terminate batch scan

        
    def nextScanName(self, scname):
        nextscname = '0000'
        if len(scname) > 0:
            scnumber = int(scname[7:11])
            nextsc_str = str(scnumber+1).zfill(4)
            nextscname = scname
            nextscname = scname.replace(scname[7:11], nextsc_str)
            
        return nextscname
        
        
    def startScan(self):
        status = 1
        scan_mode = self.scandic['scanMode']
        scan_setting = self.scandic[scan_mode]
        scanparm = []
        
        if scan_setting['pre_parm'] is not None:
            for s_ in scan_setting['pre_parm']:
                p_ = scan_setting[s_]
                scanparm += p_
        
        parm_labels = scan_setting['parm_label'] 
        scans = scan_setting['scans']

        for scan_idx, scan_ in enumerate(scans):
            t = self.getCurrentTime()
            self.logger('%s: Setting up %d/%d batch scan using %s mode.\n'%(t,
                 scan_idx+1, len(scans), scan_mode))
            next_sc = self.nextScanName(caget(self.pvs['fname_saveData']))
            self.logger('%s Initiating scan %s %s\n'%('#'*20, next_sc, '#'*20))
            self.logger('Sample info: %s\n'%(self.scandic['smpInfo']))
#             self.logger('Sample temp (K): %.3f\n'%(caget(self.pvs['temp'])))
#             logpvs = ['CryoCon1:In_3', 'CryoCon1:In_2', 'CryoCon3:In_2', 'CryoCon3:Loop_2']
#             for lpv in logpvs:
#                 self.logger('%s: %.3f\n'%(lpv, caget(self.pvs[lpv])))

            params = []
            if type(scan_) is not list:
                params = scanparm + [scan_]
            else:
                params = scanparm + scan_

            self.logger('%s: '%(self.getCurrentTime()))
            for l_, p_ in zip(parm_labels, params):
                self.logger('%s: %.3f \t'%(l_, p_))
            self.logger('\n\n')

            if scan_mode == 'angleSweep':
                self.angleSweepScanInit(parm_labels, params, next_sc, scan_idx, len(scans))
            elif scan_mode == 'fine':
                self.fineScanInit(parm_labels, params, next_sc, scan_idx, len(scans))
            elif scan_mode == 'batchXRF_fixAngle':
                self.batchXRFInit(parm_labels, params, next_sc, scan_idx, len(scans))
            else:
                status = self.coarseFineScanInit(parm_labels, params, next_sc, 
                                                 scan_idx, len(scans), scan_setting)
                
            if status == -1:
                break

        self.changeXYcombinedMode()
        if status:          
            self.logger('%s: Complete. Congratulation!\n'%(self.getCurrentTime()))
        else:
            self.logger('%s: Batch scan termiinated\n'%(self.getCurrentTime()))


        self.logfid.close()
            




