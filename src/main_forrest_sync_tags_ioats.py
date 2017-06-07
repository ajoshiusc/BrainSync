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
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt

p_dir_ref = '/big_disk/ajoshi/HCP_data/'
hemi = 'left'
ref = '100307'
TR = 2
fmri_run3 = loadmat('/deneb_disk/studyforrest/sub-02-run3\
/fmri_tnlm_0p5_reduce3_v2.mat')  # h5py.File(fname1);

dfs_ref = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.smooth.' + hemi + '.dfs'))

segl = 219
fseg1 = normdata(fmri_run3['func_'+hemi][:, 4:segl])
fseg2 = normdata(fmri_run3['func_'+hemi][:, segl:2*segl-4])

annot = pd.read_csv('/deneb_disk/studyforrest/ioats_2s_av_allchar.csv')

face_annot = sp.array(annot['face'])

tst = int(1760.0/TR)
tend = int(2620.0/TR)
face_run3 = face_annot[tst:tend]
faceseg1 = face_run3[:fseg1.shape[1]]
faceseg2 = face_run3[fseg1.shape[1]:]
faceseg1 = normdata(faceseg1[None, :]).squeeze()
faceseg2 = normdata(faceseg2[None, :]).squeeze()

fseg1 = normdata(fseg1)
fseg2 = normdata(fseg2)

#fn = sp.load('movie_corr.npz')
#rho_direct22 = fn['rho_direct22']
#ind = rho_direct22 > 0.4
#fseg1 = fseg1[ind, :]
#fseg2 = fseg2[ind, :]

fseg1_2, R = rot_sub_data(ref=fseg2, sub=fseg1)

faceseg1_2 = sp.dot(faceseg1, R.T)

faceseg1_2 = gaussian_filter(faceseg1_2, 4)

faceseg1 = gaussian_filter(faceseg1, 6)/1.4
faceseg2 = gaussian_filter(faceseg2, 6)/1.4
faceseg1 = faceseg1[:180]
faceseg2 = faceseg2[:180]
faceseg1_2 = faceseg1_2[:180]

print(sp.linalg.norm(fseg1-fseg2), sp.linalg.norm(fseg1_2-fseg2), sp.linalg.norm(fseg1_2-fseg1))
print(sp.dot(faceseg1,faceseg2)/len(faceseg2), sp.dot(faceseg1_2,faceseg2)/len(faceseg2), sp.dot(faceseg1_2,faceseg1)/len(faceseg2))

plt.plot(faceseg1, 'b')
plt.plot(faceseg2, 'r')
plt.plot(faceseg1_2, 'k')

plt.savefig('face_annotation_sync.png')
