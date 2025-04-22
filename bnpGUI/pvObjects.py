"""
Created on Wed Oct 27 10:28:14 2021

@author: graceluo

Create and get PVobjects

"""
#!/home/beams/USERBNP/.conda/envs/py36/bin/python

import epics, sys, datetime
import numpy as np
from misc import getCurrentTime
import epics.devices
from epics import caput, caget

# Eiger object is create for accessing camera related attributes
class eiger(object):
    def __init__(self, cam_pv_str, file_pv_str):
        self.pvstr = cam_pv_str
        self.cam = epics.devices.AD_Camera(cam_pv_str)
        self.fileIO = epics.devices.AD_FilePlugin(file_pv_str)
        
    def setNumTriggers(self, numTriggers):
        caput('%sNumTriggers'%self.pvstr, numTriggers)
    
    def getNumTriggers(self):
        return caget('%sNumTriggers'%self.pvstr)


class pvObject(object):
    def __init__(self, pv_str, pv_key, onchange_callback=False):
        self.pv = epics.PV(pv_str)
        self.pvname = pv_key
        self.putvalue = self.pv.value
        self.put_complete = 0
        self.motor_ready = 1
        self.time_pre = None    #datetime of PV when connected or previous value change
        self.time_delta = 0    #time difference in sec btw its value change
        if onchange_callback:
            self.pv.add_callback(self.onChanges)
        
            
    def onPutComplete(self, pvname=None,  **kws):
        sys.stdout.write('%s: Finish updating PV %s with value of %s\n'\
                          %(getCurrentTime(), self.pvname, str(self.putvalue)))
        self.put_complete = 1
        
    def onChanges(self, pvname=None, **kws):
            
        if self.time_pre is None: 
            self.time_pre = datetime.datetime.now()
        else:
            curtime = datetime.datetime.now()
            self.time_delta = (curtime-self.time_pre).seconds
            self.time_pre = curtime
            
        sys.stdout.write('%s: previous time:%s, delta time:%s\n'
                         %(getCurrentTime(), self.time_pre, self.time_delta))

        
    def put_callback(self, v = None):
        self.put_complete = 0
        if v is not None:
            self.putvalue = v
            self.pv.put(self.putvalue, callback=self.onPutComplete)

    def motorReady(self, rqspv, tolerance = 4e-2):
        rqsvalue = np.round(rqspv.value, 2)
        if abs((np.round(self.pv.value, 2) - rqsvalue)) < tolerance:
            self.motor_ready = 1
        else:
            rqspv.put(rqsvalue)
            self.motor_ready = 0
        
        
def definePVs():
     pvs = {'x_center_Rqs':'9idbTAU:SM:PX:RqsPos', 'x_center_Act':'9idbTAU:SM:PX:ActPos',
            'y_center_Rqs':'9idbTAU:SY:PY:RqsPos', 'y_center_Act':'9idbTAU:SY:PY:ActPos',
            'z_value_Rqs':'9idbTAU:SM:SZ:RqsPos', 'z_value_Act':'9idbTAU:SM:SZ:ActPos',
            'tomo_rot_Rqs':'9idbTAU:SM:CT:RqsPos', 'tomo_rot_Act':'9idbTAU:SM:CT:ActPos',
            'sm_rot_Rqs':'9idbTAU:SM:ST:RqsPos', 'sm_rot_Act':'9idbTAU:SM:ST:ActPos',
            'x_width':'9idbBNP:scan1.P1WD', 'y_width':'9idbBNP:scan2.P1WD',
            'x_step':'9idbBNP:scan1.P1SI', 'y_step':'9idbBNP:scan2.P1SI',
            'dwell':'9idbBNP:scanTran3.C', 'BDA_pos':'9idbTAU:UA:UX:RqsPos',
            'det_time':'9idbBNP:3820:ElapsedReal', '1D_time':'9idbBNP:scanTran4.F',
            'xmap_stp':'9idbXMAP:StopAll', 'netCDF_stp':'9idbXMAP:netCDF1:Capture',
            'mcs_stp':'9idbBNP:3820:StopAll', 'mcs_status':'9idbBNP:3820:Acquiring',
            'xmap_status':'9idbXMAP:Acquiring', 'netCDF_save':'9idbXMAP:netCDF1:WriteFile',
            'netCDF_status':'9idbXMAP:netCDF1:WriteFile_RBV',
            'collect_mode':'9idbXMAP:CollectMode',
            'y_motor_ready':'9idbTAU:SY:Ps:Ready', 'xztp_motor_ready':'9idbTAU:SM:Ps:Ready',
            'x_piezo_val':'9idbTAU:M7009.VAL', 'y_piezo_val':'9idbTAU:M7010.VAL',
            'scan2Record':'9idbBNP:scan2',
            
            'mono_mode': '9idb:mono_pid1.FBON', 'read_1':'9idbXMAP:scan1.R1PV', 
            'drive_1':'9idbXMAP:scan1.P1PV', 'mono_eng':'2ida2:BraggEAO.VAL',
            'dwell_step': '9idbXMAP:userTran1.P', 'xanes_eng_cen':'9idbXMAP:scan1.P1CP',
            
            'x_motorMode':'9idbTAU:SM:Ps:xMotionChoice.VAL',
            'y_motorMode':'9idbTAU:SY:Ps:yMotionChoice.VAL',
            'x_updatecenter':'9idbBNP:scan1.P1CP', 'y_updatecenter':'9idbBNP:scan2.P1CP',
            # 'x_setcenter':'9idbBNP:aoRecord11.PROC', 'y_setcenter':'9idbBNP:aoRecord12.PROC',
            'piezo_xCenter':'9idbTAU:SM:Ps:xCenter.PROC',
            'piezo_yCenter':'9idbTAU:SY:Ps:yCenter.PROC',
            'tot_lines':'9idbBNP:scan2.NPTS', 'cur_lines':'9idbBNP:scan2.CPT',
            'tot_pts_perline':'9idbBNP:scan1.NPTS',
            

            'CryoCon1:In_1':'9idbCRYO:CryoCon1:In_1:Temp.VAL',
            'CryoCon1:In_3':'9idbCRYO:CryoCon1:In_3:Temp.VAL',
            'CryoCon1:In_2':'9idbCRYO:CryoCon1:In_2:Temp.VAL',
            'CryoCon3:In_2':'9idbCRYO:CryoCon3:In_2:Temp.VAL',
            'CryoCon3:Loop_2':'9idbCRYO:CryoCon3:Loop_2:SetControl.VAL',


            'run':'9idbBNP:scan2.EXSC', 'wait':'9idbBNP:scan2.WAIT', 'wait_val':'9idbBNP:scan2.WCNT',
            'pause':'9idbBNP:scan1.PAUS', 'abort':'9idbBNP:AbortScans.PROC', 
            'msg1d':'9idbBNP:scan1.SMSG',
            'fname_saveData':'9idbBNP:saveData_fileName',
            'filesys':'9idbBNP:saveData_fileSystem',
            'subdir':'9idbBNP:saveData_subDir',
            'nextsc':'9idbBNP:saveData_scanNumber',
            'basename':'9idbBNP:saveData_baseName',
            
            }
     return pvs
 
def scan2RecordDetectorTrigerPVs():
    pvs = {'scan1':'9idbBNP:scan1.EXSC',
           'eigerAcquire':'2iddEGR:cam1:Acquire',
           'eigerFileCapture':'2iddEGR:HDF1:Capture'}
    return pvs
    
def getEiger():
    # create Eiger cam record
    e = eiger('2iddEGR:cam1:', '2iddEGR:HDF1:')
    return e
    
#    pvs = {'test1':'2idbleps:userTran2.CMTA', 'test2':'2idbleps:userTran2.CMTB', 'test3':'2idbleps:userTran2.CMTC',
#           'test4':'2idbleps:userTran2.CMTD', 'test5':'2idbleps:userTran2.CMTE', 'test6':'2idbleps:userTran2.CMTF',
#           'test7':'2idbleps:userTran2.CMTG'}


def getPVobj():
    pvObjs = {}
    pvs = definePVs()
    for k, v in pvs.items():
        if 'Record' not in k:
            pv_obj = pvObject(v, k, onchange_callback=True if '9idbBNP:scan2.CPT'==v else False)
            pvObjs.update({k: pv_obj})
        else:
            pvObjs.update({k:epics.devices.Scan(v)})
    return pvObjs
