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
from fmri_methods_sipi import rot_sub_data, reorder_labels
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
hemi = 'right'
ref = '196750'
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
dfs_hemi = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.' + hemi + '.dfs'))
dfs_hemi_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.very_smooth.' + hemi + '.dfs'))

surf1 = dfs_hemi_sm

lst = os.listdir('/big_disk/ajoshi/HCP5')
rho1 = 0
rho1rot = 0
rho2 = 0
rho2rot = 0
# lst = [lst[0]]
diffbefore = 0
diffafter = 0

sub = lst[0]

vrest1 = scipy.io.loadmat('/big_disk/ajoshi/coding_ground/brainsync/data/\
Cleveland/subject1/fmri_tnlm_5_reduce3_v2.mat')  # h5py.File(fname1);


data = vrest1['func_' + hemi + '']
indx = sp.isnan(data)
data[indx] = 0


vrest = data
m = np.mean(vrest, 1)
vrest = vrest - m[:, None]
''' This is added since length of fcon1000 sequences is 225'''
#vrest = vrest[:, :225]

s = np.std(vrest, 1)+1e-116
vrest1 = vrest/s[:, None]

rho1 = 0
rho1rot = 0
diffafter = 0
diffbefore = 0

a = sp.load('/big_disk/ajoshi/coding_ground/brainsync/data/\
fcon1000_null_all_' + hemi + '.npz')
rho_null = sp.mean(a['rho_null'], axis=0)

lst = glob.glob('/big_disk/ajoshi/fcon_1000/Beijing/sub*')
nsub = 0
rho_all = sp.zeros((vrest1.shape[0], 0))



for sub in lst:
    if not os.path.exists(sub + '/fmri_tnlm_5_reduce3_v2.mat'):
        continue

    vrest2 = scipy.io.loadmat(sub + '/fmri_tnlm_5_reduce3_v2.mat')        
    data = vrest2['func_' + hemi + '']
    indx = sp.isnan(data)
    data[indx] = 0
    vrest = data
    vrest = vrest[:, :vrest1.shape[1]]
    m = np.mean(vrest, 1)
    vrest = vrest - m[:, None]
    s = np.std(vrest, 1)+1e-116
    vrest2 = vrest/s[:, None]

    rho1 += sp.sum(vrest1*vrest2, axis=1)/vrest1.shape[1]
    diffbefore += vrest1 - vrest2

    vrest2, Rot, _ = rot_sub_data(ref=vrest1, sub=vrest2,
                                  area_weight=sp.sqrt(surf_weight))
    t = sp.sum(vrest1*vrest2, axis=1)/vrest1.shape[1]

    rho_all = sp.append(rho_all, t[:, None], axis=1)
    rho1rot += t

    diffafter += vrest1 - vrest2
    nsub += 1
    print(sub)

rho1rot /= nsub


#
# import seaborn as sns
#
# sns.distplot(rho_all[120,:])
# sns.distplot(rho_null[120,:])

#
# pval=sp.zeros((rho_all.shape[0],1))
#
# for jj in range(rho_all.shape[0]):
#     _, pval[jj] = sp.stats.mannwhitneyu(rho_null[jj,:], rho_all[jj,:]) #,
#                                         # alternative='greater')
#     print jj
# sns.distplot(pval)
dfs_hemi_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.smooth.' + hemi + '.dfs'))

dfs_hemi_sm.attributes = sp.squeeze(rho_all.mean(axis=1))
dfs_hemi_sm = patch_color_attrib(dfs_hemi_sm, clim=[0, 1])
view_patch_vtk(dfs_hemi_sm, azimuth=90, elevation=180, roll=90,
               outfile='rest_rot1_fcon1000_subject1_' + hemi + '.png', show=1)
view_patch_vtk(dfs_hemi_sm, azimuth=-90, elevation=180, roll=-90,
               outfile='rest_rot2_fcon1000_subject1_' + hemi + '.png', show=1)

rho_null = rho_null.T
rho_all1 = sp.mean(rho_all, axis=1)[:, None]
pval = 0*rho_all1

for jj in range(rho_all1.shape[1]):
    pval[:, jj] = sp.sum(rho_null > rho_all1[:, jj][:, None], axis=1)

pval1 = 1-pval.mean(axis=1)/rho_null.shape[1]
a = multipletests(pvals=pval1, alpha=0.05, method='fdr_bh')
dfs_hemi_sm.attributes = pval1 #a[1]+1 # sp.amax((pval1, a[1]),axis=0)
# sns.distplot(pval1)

# from rpy2.robjects.packages import importr
# from rpy2.robjects.vectors import FloatVector
#
# stats = importr('stats')
#
# p_adjust = stats.p_adjust(FloatVector(pval1), method = 'BH')
# dfs_hemi_sm.attributes= sp.array(p_adjust)

dfs_hemi_sm = patch_color_attrib(dfs_hemi_sm, clim=[0, 0.05])
view_patch_vtk(dfs_hemi_sm, azimuth=90, elevation=180, roll=90,
               outfile='rest_after_rot1_fcon1000_subject1_' + hemi + '_fdr.png', show=1)
view_patch_vtk(dfs_hemi_sm, azimuth=-90, elevation=180, roll=-90,
               outfile='rest_after_rot2_fcon1000_subject1_' + hemi + '_fdr.png', show=1)
