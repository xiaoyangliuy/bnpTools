"""
Created on Wed Jul 28 14:32:26 2021

@author: graceluo

Function to setup scan, interact with PV mostly
"""

from pvObjects import getPVobj
from misc import getCurrentTime
import os, time, sys
import numpy as np

class pvComm():
    def __init__(self, userdir = None, log = 'log.txt'):
        self.pvs = getPVobj()
        if userdir is None:
            self.userdir = self.getDir()
        else:
            self.userdir = userdir
        self.logfilepath = os.path.join(self.userdir, log)
        self.logfid = open(self.logfilepath, 'a')
            
    def logger(self, msg):
        sys.stdout.write(msg)
        sys.stdout.flush()
        if self.logfid.closed:
            self.logfid = open(self.logfilepath, 'a')
        self.logfid.write(msg)
        self.logfid.flush()
    
    def getDir(self):
        fs = self.pvs['filesys'].pv.value
        fs = fs.replace('//micdata/data1', '/mnt/micdata1')
        return os.path.join(fs, self.pvs['subdir'].pv.value.replace('mda', ''))
    
    def getBDAx(self):
        return np.round(self.pvs['BDA_pos'].pv.value, 2)
    
    def getSMAngle(self):
        return np.round(self.pvs['sm_rot_Act'].pv.value, 2)
    
    def getTomoAngle(self):
        return np.round(self.pvs['tomo_rot_Act'].pv.value, 2)
    
    def scanPause(self):
        self.pvs['wait'].put_callback(1)
    
    def scanResume(self):
        self.pvs['wait'].put_callback(0)
        
    def scanAbort(self):
        self.pvs['abort'].put_callback(1)
    
    def changeTomoRotate(self, theta):
        curr_angle = np.round(self.pvs['tomo_rot_Act'].pv.value, 2)
        t = getCurrentTime()
        self.logger('%s; Changing tomo rotation angle from to %.2f to %.2f\n'%(t, curr_angle, theta))
        self.pvs['tomo_rot_Act'].put_callback(theta)
    
    def changeSMRotate(self, theta):
        curr_angle = np.round(self.pvs['sm_rot_Act'].curr_value, 2)
        t = getCurrentTime()
        self.logger('%s; Changing sample rotation angle from to %.2f to %.2f\n'%(t, curr_angle, theta))
        self.pvs['sm_rot_Act'].put_callback(theta)
    
    def blockBeamBDA(self, BDA):
        bda_pos = BDA - 500
        t = getCurrentTime()
        self.logger('%s: Move BDA to block position at: %.3f\n'%(t, bda_pos))
        self.pvs['BDA_pos'].put_callback(bda_pos)
        
    def openBeamBDA(self, BDA):
        self.logger('%s: Move BDA to open position at: %.3f\n'%(getCurrentTime(), BDA))
        self.pvs['BDA_pos'].put_callback(BDA)
    
    def changeXYcombinedMode(self):        
        self.logger('%s; Changing XY scan mode to combined motion\n'%(getCurrentTime()))
        self.pvs['x_motorMode'].pv.put(0)
        self.pvs['y_motorMode'].pv.put(0)
        
    def changeXtoCombinedMode(self):        
        self.logger('%s; Changing XY scan mode to combined motion\n'%(getCurrentTime()))
        self.pvs['x_motorMode'].pv.put(0) 
        
    def changeXtoPiezolMode(self):
        self.logger('%s: Changing X scan mode to Piezo only\n'%(getCurrentTime()))
        self.pvs['x_motorMode'].pv.put(2)

    def setXYcenter(self):
        self.logger('%s: Update the current position as the center of'\
                    'the scan.\n'%(getCurrentTime()))
        self.pvs['x_setcenter'].pv.put(1)
        self.pvs['y_setcenter'].pv.put(1)
        
    def centerPiezoXY(self):
        self.logger('%s: Centering piezoX and piezoY.\n'%(getCurrentTime()))
        self.pvs['piezo_xCenter'].pv.put(1)
        time.sleep(1)
        self.pvs['piezo_yCenter'].pv.put(1)
        time.sleep(1)
    
    def assignPosValToPVs(self, pvstr, pvval):
        for s_, v_ in zip(pvstr, pvval):
            self.pvs[s_].pv.put(v_)
            self.logger('%s: Change %s to %.3f\n' % (getCurrentTime(), s_, v_))
            
    def assignSinglePV(self, pvstr, pvval):
        self.pvs[pvstr].pv.put(pvval)
        self.logger('%s: Change %s to %.3f\n' % (getCurrentTime(), pvstr, pvval))
            
    def writeScanInit(self, mode, smpinfo, scandic):
        next_sc = self.nextScanName()
        self.logger('%s Initiating scan %s %s\n'%('#'*20, next_sc, '#'*20))
        self.logger('Sample info: %s\n'% smpinfo)
        self.logger('%s: Setting up scan using %s mode.\n'%(getCurrentTime(), mode))
        self.logger('%s: %s'%(getCurrentTime(), scandic))
        self.logger('\n\n')
        
    def motorReady(self, l, mt):
        self.logger('%s: Checking whether motors are ready.\n'%(getCurrentTime()))
        actpv = self.pvs['%s_Act'%l].pv
        rqspv = self.pvs['%s_Rqs'%l].pv
        self.pvs['%s_Act'%l].motorReady(rqspv, mt)
        
        if self.pvs['%s_Act'%l].motor_ready:
            self.logger('%s: %s motor is in position with value'\
                            '%.2f\n'%(getCurrentTime(), l, actpv.value))
            return 1
        else:
            self.logger('%s: %s motor not in position, current: %.2f,'\
                            ' request: %.2f\n'%(getCurrentTime(), l, actpv.value, rqspv.value))
            return 0
    
#    def motorReady(self, label, mtolerance):
#        self.logger('%s: Checking whether motors are ready.\n'%(getCurrentTime()))
#        rules = [0] * len(label)
#        for i, (l, t) in enumerate(zip(label, mtolerance)):
#            actpv = self.pvs['%s_Act'%l].pv
#            rqspv = self.pvs['%s_Rqs'%l].pv
#            self.pvs['%s_Act'%l].motorReady(rqspv, t)
#            if self.pvs['%s_Act'%l].motor_ready:
#                self.logger('%s: %s motor is in position with value'\
#                                '%.2f\n'%(getCurrentTime(), l, actpv.value))
#                rules[i] = 1
#            else:
#                self.logger('%s: %s motor not in position, current: %.2f,'\
#                                ' request: %.2f\n'%(getCurrentTime(), l, actpv.value, rqspv.value))
#            time.sleep(1)
#                
#        if all(rules):
#            self.logger('%s: Motors ready\n'%getCurrentTime())
#            return 1
#        else:
#            return 0         
    
    def nextScanName(self):
        return '%s%s.mda'%(self.pvs['basename'].pv.value, 
                           str(self.pvs['nextsc'].pv.value).zfill(4))
    
    def getXYZcenter(self):
        return [np.round(self.pvs[i].pv.value, 2) for i in ['x_center_Act', \
                'y_center_Act', 'z_value_Act']]
        

    

    

