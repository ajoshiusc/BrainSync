# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 05:39:43 2016

@author: ajoshi
"""
import sys
import scipy.io
import scipy as sp
import os
import numpy as np
import nibabel as nib
from dfsio import readdfs, writedfs
from surfproc import view_patch, view_patch_vtk, smooth_surf_function,\
                     face_v_conn, patch_color_attrib
from fmri_methods_sipi import rot_sub_data, reorder_labels, normdata
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import glob
from statsmodels.stats.multitest import multipletests

p_dir = '/big_disk/ajoshi/HCP_data'
p_dir_ref = '/big_disk/ajoshi/HCP_data/'
lst = os.listdir(p_dir)
r_factor = 3
ref_dir = os.path.join(p_dir_ref, 'reference')
nClusters = 3

ref = '100307'
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
dfs_left = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.left.dfs'))
dfs_left_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.very_smooth.left.dfs'))

surf1 = dfs_left_sm

# smooth_surf_function(dfs_left_sm, Wt*Ar*0.1, a1=0, a2=1)

# sub = '110411'
# p_dir = '/home/ajoshi/data/HCP_data'
lst = os.listdir('/big_disk/ajoshi/HCP5')
rho1 = 0
rho1rot = 0
rho2 = 0
rho2rot = 0
# lst = [lst[0]]
diffbefore = 0
diffafter = 0

sub = lst[0]

vrest1 = scipy.io.loadmat('/big_disk/ajoshi/fcon_1000/Beijing/sub74386/fmrit_reduce3_v2.mat')  # h5py.File(fname1);
data = vrest1['fmri_left']
vrest1 = normdata(data)

vrest2 = scipy.io.loadmat('/big_disk/ajoshi/fcon_1000/Beijing/sub56136/fmrit_reduce3_v2.mat')  # h5py.File(fname1);
data = vrest2['fmri_left']
vrest2d = normdata(data)
synced={}
synced['synced_left'],_ = rot_sub_data(ref=vrest1, sub=vrest2d)
synced['left1']=vrest1
synced['left2']=vrest2d

vrest1 = scipy.io.loadmat('/big_disk/ajoshi/fcon_1000/Beijing/sub74386/fmrit_reduce3_v2.mat')  # h5py.File(fname1);
data = vrest1['fmri_right']
vrest1 = normdata(data)

vrest2 = scipy.io.loadmat('/big_disk/ajoshi/fcon_1000/Beijing/sub56136/fmrit_reduce3_v2.mat')  # h5py.File(fname1);
data = vrest2['fmri_right']
vrest2d = normdata(data)
synced['synced_right'],_ = rot_sub_data(ref=vrest1, sub=vrest2d)
synced['right1']=vrest1
synced['right2']=vrest2d



sp.io.savemat('out_synced.mat',synced)