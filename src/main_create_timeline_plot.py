#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 18:53:24 2017

@author: ajoshi
"""

import scipy.io
import scipy as sp
import os
import numpy as np
from dfsio import readdfs
from surfproc import view_patch_vtk, patch_color_attrib
#from fmri_methods_sipi import rot_sub_data
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

s = sp.zeros(284)
# LR Motor
lh = sp.int16([71.518/.72, (71.518+12)/.72, 162.015/0.72,
              (162.015 + 12)/0.72])
rh = sp.int16([11.009/.72, (11.009+12)/.72, 131.894/0.72,
              (131.894 + 12)/0.72])
lf = sp.int16([26.136/.72, (26.136+12)/.72, 116.766/0.72,
              (116.766 + 12)/0.72])
rf = sp.int16([56.391/.72, (56.391+12)/.72, 177.142/0.72,
              (177.142 + 12)/0.72])
t = sp.int16([41.264/.72, (41.264+12)/.72, 101.639/0.72,
             (101.639 + 12)/0.72])

s[lh[0]:lh[1]] = 1
s[lh[2]:lh[3]] = 1
s[rh[0]:rh[1]] = -1
s[rh[2]:rh[3]] = -1
s[lf[0]:lf[1]] = 2
s[lf[2]:lf[3]] = 2
s[rf[0]:rf[1]] = -2
s[rf[2]:rf[3]] = -2
s[t[0]:t[1]] = 3
s[t[2]:t[3]] = 3
#
## RL Motor
#lh = 284+sp.int16([11.009/.72, (11.009+12)/.72, 116.633/0.72,
#                  (116.633 + 12)/0.72])
#rh = 284+sp.int16([86.512/.72, (86.512+12)/.72, 162.014/0.72,
#                  (162.014 + 12)/0.72])
#lf = 284+sp.int16([71.384/.72, (71.384+12)/.72, 177.141/0.72,
#                  (177.141 + 12)/0.72])
#rf = 284+sp.int16([26.136/.72, (26.136+12)/.72, 146.887/0.72,
#                  (146.887 + 12)/0.72])
#t = 284+sp.int16([56.257/.72, (56.257+12)/.72, 131.76/0.72,
#                 (131.76+12)/0.72])
#
#s[lh[0]:lh[1]] = 1
#s[lh[2]:lh[3]] = 1
#s[rh[0]:rh[1]] = -1
#s[rh[2]:rh[3]] = -1
#s[lf[0]:lf[1]] = 2
#s[lf[2]:lf[3]] = 2
#s[rf[0]:rf[1]] = -2
#s[rf[2]:rf[3]] = -2
#s[t[0]:t[1]] = 3
#s[t[2]:t[3]] = 3

so = s.copy()

from pylab import rcParams
rcParams['figure.figsize'] = 20, 2
for ind in sp.arange(284):
#    plt.figure(figsize=(1200,20))
    s = so.copy()
    s[ind] = 4
    plt.plot(s)
    tname = 'timing_%d.png' % ind
    plt.savefig(tname)
    plt.close()

motor_timing = so
sp.savez('timing.npz', motor_timing)

