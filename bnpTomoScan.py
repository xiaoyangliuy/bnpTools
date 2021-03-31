'''
This is the scan class for tomo data collection at BNP. 

TO DO: add temp PV and output to logbook; remove not important lines in logbook
'''

import time, os, sys
import tqdm
from epics import caput, caget
from time import gmtime, strftime
import numpy as np
import imgProcessing

class bnpTomoScan():
    
    def __init__(self, scandic):
        self.logfid = open(os.path.join(scandic['userdir'], scandic['logfile']), 'a+')
        self.scandic = scandic
        self.pvs = self.definePVs()
        
    def logger(self, msg):
        sys.stdout.write(msg)
        sys.stdout.flush()
        self.logfid.write(msg)
        self.logfid.flush()
        
    def getCurrentTime(self):
        return strftime('%Y-%m-%d %H:%M:%S', gmtime())
    
    def definePVs(self):
        self.logger('%s: Associate motors with PVs\n'%(self.getCurrentTime()))
        pvs = {'x_center_Rqs':'9idbTAU:SM:PX:RqsPos', 'x_center_Act':'9idbTAU:SM:PX:ActPos',
               'y_center_Rqs':'9idbTAU:SY:PY:RqsPos', 'y_center_Act':'9idbTAU:SY:PY:ActPos',
               'z_value_Rqs':'9idbTAU:SM:SZ:RqsPos', 'z_value_Act':'9idbTAU:SM:SZ:ActPos',
               'tomo_rot_Rqs':'9idbTAU:SM:CT:RqsPos', 'tomo_rot_Act':'9idbTAU:SM:CT:ActPos',
               'sm_rot_Rqs':'9idbTAU:SM:ST:RqsPos', 'sm_rot_Act':'9idbTAU:SM:ST:ActPos',
               'x_width':'9idbBNP:scan1.P1WD', 'y_width':'9idbBNP:scan2.P1WD',
               'x_step':'9idbBNP:scan1.P1SI', 'y_step':'9idbBNP:scan2.P1SI',
               'dwell':'9idbBNP:scanTran3.C', 'BDA_pos':'9idbTAU:UA:UX:RqsPos',
               
               'x_motorMode':'9idbTAU:SM:Ps:xMotionChoice.VAL', 'y_motorMode':'9idbTAU:SY:Ps:yMotionChoice.VAL',
               'x_setcenter':'9idbBNP:aoRecord11.PROC', 'y_setcenter':'9idbBNP:aoRecord12.PROC',
               'piezo_xCenter':'9idbTAU:SM:Ps:xCenter.PROC', 'piezo_yCenter':'9idbTAU:SY:Ps:yCenter.PROC',
               'tot_lines':'9idbBNP:scan2.NPTS', 'cur_lines':'9idbBNP:scan2.CPT',
               'temp':'9idbCRYO:CryoCon1:In_3:Temp.VAL',
               
               'run':'9idbBNP:scan2.EXSC',
               'fname_saveData':'9idbBNP:saveData_fileName'
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
        self.logger('%s: Update the current position as the center of the scan.\n'%(self.getCurrentTime()))
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
            sys.stdout.write('%s: Motors not in position\n'%(self.getCurrentTime()))
            for i, l_ in enumerate(label):
                caput(self.pvs['%s_Rqs'%(l_)], caget(self.pvs['%s_Rqs'%(l_)]))
                time.sleep(1)
                pos_diff = abs(caget(self.pvs['%s_Rqs'%(l_)]) - caget(self.pvs['%s_Act'%(l_)]))
                rules[i] = (pos_diff < mtolerance[i])
                if rules[i]:
                    self.logger('%s: %s motor is in position with value %.2f um\n'%(self.getCurrentTime(), l_, caget(self.pvs['%s_Rqs'%(l_)])))
            
            if all(rules):
                ready = 1
                    
        self.logger('%s: Motors Ready \n'%(self.getCurrentTime()))

    def execScan(self, scname, scidx, n_scns):
        self.setXYcenter()
        self.changeXtoPiezolMode()
        self.centerPiezoXY()
        self.openBeamBDA()
        
        # Execute Scan
        caput(self.pvs['run'], 1)
        self.logger('%s: Scanning \n'%(self.getCurrentTime()))
        nlines = caget(self.pvs['tot_lines'])
        while caget(self.pvs['run']):
            time.sleep(10)
            cline = caget(self.pvs['cur_lines'])
            sys.stdout.write('Scanning %s (batch %d/%d): line %d/%d is done\n'%(scname, scidx+1, n_scns, cline, nlines))
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
    
    def coarseScanInit(self, params_label, params, scname, scidx, n_scns):
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
        
#         if self.scandic['coarse_fine']['angle0_afScan']:
#             self.logger('%s: Putting Tomo angle rotation back to 0 deg after scan\n'%(self.getCurrentTime()))
#             self.changeTomoRotate(0) # initialize angle at 0 deg
        
    def nextScanName(self, scname):
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
            self.logger('Sample temp (K): %.3f\n'%(caget(self.pvs['temp'])))
            
            params = []
            if type(scan_) is not list:
                params = scanparm + [scan_]
            else:
                params = scanparm + scan_
            
            self.logger('%s: '%(self.getCurrentTime()))
            for l_, p_ in zip(parm_labels, params):
                self.logger('%s: %.3f \t'%(l_, p_))
            self.logger('\n\n')
                
            if scan_mode == 'coarse':
                self.coarseScanInit(parm_labels, params, next_sc, scan_idx, len(scans))
            elif scan_mode == 'fine':
                self.fineScanInit(parm_labels, params, next_sc, scan_idx, len(scans))
            elif scan_mode == 'batchXRF_fixAngle':
                self.batchXRFInit(parm_labels, params, next_sc, scan_idx, len(scans))
            else:
                self.coarseScanInit(parm_labels, params, next_sc, scan_idx, len(scans))
                # If find_bbox is on, perform image processing on the coarse scan to get coordinates
                if scan_setting['find_bbox']:
                    cscan_path = os.path.join(self.scandic['userdir'], 'img.dat/%s.h5'%(next_sc))
                    while not os.path.exists(cscan_path):
                        time.sleep(1)
                    
                    time_diff = 0
                    while (time_diff < 10):
                        time.sleep(1)
                        file_mod_time = os.stat(cscan_path).st_mtime
                        time_diff = int(time.time() - file_mod_time)
                        sys.stdout.write('Waiting for coarse scan file %s to be ready, file modified time: %d, time difference: %d \n'%(next_sc, file_mod_time, time_diff))
                        
                        
                    img_path= os.path.join(self.scandic['userdir'],'imgProg')
                    if not os.path.exists(img_path):
                        os.makedirs(img_path)
                        
                    figpath = os.path.join(img_path, 'bbox_%s.png'%(next_sc))
                    new_x, new_y, new_w, new_h = imgProcessing.getROIcoordinate(cscan_path, 
                                                                                self.scandic['elm'], figpath = figpath)
                    f_scanparm = scan_setting['fine_pts_area']
                    arealim = scan_setting['area_lim']
                    if (new_w < arealim * f_scanparm[0]) | (new_y < arealim * f_scanparm[1]):
                        self.logger('Extract ROI corr from %s:\n %.2f(x-center),\n %.2f(y-center),\n'%(next_sc, new_x, new_y))
                        self.logger('%.2f(width),\n %.2f(height)'%(new_w, new_h))
                        curr_smz = caget(self.pvs['z_value_Rqs'])
                        params = []
                        params = f_scanparm + [scan_] + [new_x, new_y, curr_smz]
                        flabels = ['x_width', 'y_width', 'x_step', 'y_step', 'dwell', 'sm_rot_Rqs',
                                      'x_center_Rqs','y_center_Rqs', 'z_value_Rqs']

                        next_sc = self.nextScanName(caget(self.pvs['fname_saveData']))
                        self.logger('%s Initiating fine scan %s %s\n'%('#'*20, next_sc, '#'*20))
                        self.fineScanInit(flabels, params, next_sc, scan_idx, len(scans))
                        
                    else:
                        self.logger('ROI area from coarse scan is bigger than the set scanning area defined in fine_pts_area\n')
                        self.logger('%s: Batch scan terminated, couldnt find proper ROI area\n'%(self.getCurrentTime()))
                        status = -1  # Flag to terminate batch scan
                        
            if status == -1:
                break
                
        self.changeXYcombinedMode()
       
        if status:          
            self.logger('%s: Complete. Congratulation!\n'%(self.getCurrentTime()))
        else:
            self.logger('%s: Batch scan termiinated\n'%(self.getCurrentTime()))
        
        
        self.logfid.close()
#         input('Press Enter to Exit')
            




