import json
import numpy as np

def loadAlignedDataFromJson(fpath, elms = None):
    with open(fpath) as json_file:
        data = json.load(json_file)

    if elms is None:
        print(data.keys())
        elms = list(data.keys())[2:]
    angles = np.array(data['angles'], dtype = 'float')
    a, h, w = np.array(data[elms[0]]).shape
    elmdata = np.zeros((len(elms), a, h, w))
    
    for i, e in enumerate(elms):
        elmdata[i,...] = np.array(data[e])
    elmdata = np.moveaxis(elmdata, 1, 2)
    p_angles = angles * np.pi / 180
    
    return {'elms':elms, 'elmdata':elmdata, 'angles_deg':angles, 'angle_rad':p_angles, 'stepsize':data['stepsize']}