# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 05:39:43 2016

@author: ajoshi
"""
from fmri_methods_sipi import normdata, rot_sub_data
from scipy.io import loadmat
import scipy as sp
from dfsio import readdfs
import os.path
from surfproc import view_patch_vtk, patch_color_attrib
import pandas as pd
p_dir_ref = '/big_disk/ajoshi/HCP_data/'
hemi = 'left'
ref = '100307'

fmri_sub1 = loadmat('/deneb_disk/studyforrest/sub-01-run1\
/fmri_tnlm_5_reduce3_v2.mat')  # h5py.File(fname1);

dfs_ref = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.smooth.' + hemi + '.dfs'))

segl = 225
sub1seg1 = normdata(fmri_sub1['func_'+hemi][:, :segl])
sub1seg2 = normdata(fmri_sub1['func_'+hemi][:, segl:2*segl])

v = pd.read_csv('/home/ajoshi/Downloads/emotions_av_1s_thr50.tsv', sep='\t')

annot = []
stseg = []
edseg = []

for ind in range(v.shape[0]):
    annot.append(v.values[ind, 2])
    stseg.append(v.values[ind, 0])
    edseg.append(v.values[ind, 1])

print annot
print sp.array(stseg)
print sp.array(edseg)


arousal = sp.zeros(len(stseg))
startTime = sp.zeros(len(stseg))
endTime = sp.zeros(len(stseg))

for ind1 in range(len(stseg)):
    ind = annot[ind1].find('arousal=')
    arousal[ind1] = float(annot[ind1][ind+8:ind+13])
    startTime[ind1] = stseg[ind1]
    endTime[ind1] = edseg[ind1]

# Form annotation time series

func_left = fmri_sub1['func_left']
arousal_annot = sp.zeros(func_left.shape[1])

TR = 2
for ind1 in range(len(stseg)):
    i = sp.arange(sp.around(startTime[ind1]/TR), sp.around(endTime[ind1]/TR))
    i = sp.int16(i)
    arousal_annot[i[i < func_left.shape[1]]] = arousal[ind1]
