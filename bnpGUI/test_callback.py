#!/home/beams/USERBNP/.conda/envs/py36/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 12:58:21 2022

Test PV call back function: onChange

@author: graceluo
"""

from pvObjects import getPVobj
import time

a = getPVobj()
for i in range(100):
    time.sleep(1)
    