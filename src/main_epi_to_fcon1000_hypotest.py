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

ref = '100307'
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
dfs_left = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.left.dfs'))
dfs_left_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.very_smooth.left.dfs'))

surf1 = dfs_left_sm
X = surf1.vertices[:, 0]
Y = surf1.vertices[:, 1]
Z = surf1.vertices[:, 2]
NumTri = surf1.faces.shape[0]
#    NumVertx = X.shape[0]
vertx_1 = surf1.faces[:, 0]
vertx_2 = surf1.faces[:, 1]
vertx_3 = surf1.faces[:, 2]
V1 = np.column_stack((X[vertx_1], Y[vertx_1], Z[vertx_1]))
V2 = np.column_stack((X[vertx_2], Y[vertx_2], Z[vertx_2]))
V3 = np.column_stack((X[vertx_3], Y[vertx_3], Z[vertx_3]))
x1 = np.zeros((NumTri))
y1 = np.zeros((NumTri))
v2_v1temp = V2-V1
x2 = np.linalg.norm(v2_v1temp, axis=1)
y2 = np.zeros((NumTri))
x3 = np.einsum('ij,ij->i', (V3-V1),
               (v2_v1temp/np.column_stack((x2, x2, x2))))
mynorm = np.cross((v2_v1temp), V3-V1, axis=1)
yunit = np.cross(mynorm, v2_v1temp, axis=1)
y3 = np.einsum('ij,ij->i', yunit, (V3-V1))/np.linalg.norm(yunit, axis=1)
sqrt_DT = (np.abs((x1*y2 - y1*x2)+(x2*y3 - y2*x3)+(x3*y1 - y3*x1)))
Ar = 0.5*(np.abs((x1*y2 - y1*x2)+(x2*y3 - y2*x3)+(x3*y1 - y3*x1)))

TC = face_v_conn(surf1)
Wt = (1.0/3.0)*(TC)
# Wt = sp.sparse.spdiags(Wt*Ar, (0), NumTri, NumTri)
surf_weight = Wt*Ar
surf1.attributes = surf_weight
surf_weight = surf_weight[:, None]
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

vrest1 = scipy.io.loadmat('/big_disk/ajoshi/coding_ground/epilepsy/\
Cleveland/subject1/fmri_tnlm_5_reduce3_v2.mat')  # h5py.File(fname1);


data = vrest1['func_left']
indx = sp.isnan(data)
data[indx] = 0


vrest = data
m = np.mean(vrest, 1)
vrest = vrest - m[:, None]
s = np.std(vrest, 1)+1e-116
vrest1 = vrest/s[:, None]

rho1 = 0
rho1rot = 0
diffafter = 0
diffbefore = 0

a = sp.load('fcon1000_null_all_left.npz')
rho_null = sp.mean(a['rho_null'], axis=0)

lst = glob.glob('/big_disk/ajoshi/fcon_1000/Beijing/sub*')
nsub = 0
rho_all = sp.zeros((vrest1.shape[0], 0))

#sub = lst[4]
#vrest2 = scipy.io.loadmat(sub + '/fmri_tnlm_5_reduce3_v2.mat')  # h5py.File(fname1);
#data = vrest2['func_left']
#indx = sp.isnan(data)
#data[indx] = 0
#vrest = data
#vrest = vrest[:, :vrest1.shape[1]]
#m = np.mean(vrest, 1)
#vrest = vrest - m[:, None]
#s = np.std(vrest, 1)+1e-116
#vrest1 = vrest/s[:, None]

for sub in lst:
    if not os.path.exists(sub + '/fmri_tnlm_5_reduce3_v2.mat'):
        continue

    vrest2 = scipy.io.loadmat(sub + '/fmri_tnlm_5_reduce3_v2.mat')
    data = vrest2['func_left']
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

    vrest2, Rot = rot_sub_data(ref=vrest1, sub=vrest2,
                               area_weight=sp.sqrt(surf_weight))
    t = sp.sum(vrest1*vrest2, axis=1)/vrest1.shape[1]

    rho_all = sp.append(rho_all, t[:, None], axis=1)
    rho1rot += t

    diffafter += vrest1 - vrest2
    nsub += 1
    print sub

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
#
dfs_left_sm.attributes = sp.squeeze(rho_all.mean(axis=1))
dfs_left_sm = patch_color_attrib(dfs_left_sm, clim=[0, 1])
view_patch_vtk(dfs_left_sm, azimuth=90, elevation=180, roll=90,
               outfile='rest_rot1_fcon1000_subject1_left.png', show=1)
view_patch_vtk(dfs_left_sm, azimuth=-90, elevation=180, roll=-90,
               outfile='rest_rot2_fcon1000_subject1_left.png', show=1)

rho_null = rho_null.T
rho_all1 = sp.mean(rho_all, axis=1)[:, None]
pval = 0*rho_all1

for jj in range(rho_all1.shape[1]):
    pval[:, jj] = sp.sum(rho_null > rho_all1[:, jj][:, None], axis=1)

pval1 = 1-pval.mean(axis=1)/rho_null.shape[1]
a = multipletests(pvals=pval1, alpha=0.05, method='fdr_bh')
dfs_left_sm.attributes = a[1] # sp.amax((pval1, a[1]),axis=0)
# sns.distplot(pval1)

# from rpy2.robjects.packages import importr
# from rpy2.robjects.vectors import FloatVector
#
# stats = importr('stats')
#
# p_adjust = stats.p_adjust(FloatVector(pval1), method = 'BH')
# dfs_left_sm.attributes= sp.array(p_adjust)

dfs_left_sm = patch_color_attrib(dfs_left_sm, clim=[0, 0.05])
view_patch_vtk(dfs_left_sm, azimuth=90, elevation=180, roll=90,
               outfile='rest_after_rot1_fcon1000_subject1_left_fdr.png', show=1)
view_patch_vtk(dfs_left_sm, azimuth=-90, elevation=180, roll=-90,
               outfile='rest_after_rot2_fcon1000_subject1_left_fdr.png', show=1)

