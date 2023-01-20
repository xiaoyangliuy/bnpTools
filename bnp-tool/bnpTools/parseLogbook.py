#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import collections


def parse_logbook(logpath):
    sc_dic = collections.defaultdict(list)
    with open(logpath, 'r') as f:
        for l in f:
            if 'Initiating scan' in l:
                smpInfo = next(f)
                smpInfo = smpInfo.replace('\n',' ')
                parmInfo = next(f)
                parmInfo = parmInfo.split('\t')
                parmInfo[0] = parmInfo[0][21:]

                sc = l[l.index('_')+1:l.index('.')]
                smpInfo = smpInfo[smpInfo.index(':')+1:]
                sc_dic['scan_number'].append(sc)
                sc_dic['smp_info'].append(smpInfo)
                for p_ in parmInfo:
                    if p_ != '\n':
                        t = p_.split(':')
                        try:
                            sc_dic[t[0]].append(float(t[1]))
                        except:
                            print(t)
    return pd.DataFrame(sc_dic)

