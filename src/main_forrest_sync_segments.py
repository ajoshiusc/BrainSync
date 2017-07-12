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

p_dir_ref = '/big_disk/ajoshi/HCP_data/'
hemi = 'left'
ref = '100307'

fmri_sub1 = loadmat('/deneb_disk/studyforrest/sub-01-run1\
/fmri_tnlm_5_reduce3_v2.mat')  # h5py.File(fname1);

fmri_sub2 = loadmat('/deneb_disk/studyforrest/sub-02-run1\
/fmri_tnlm_5_reduce3_v2.mat')  # h5py.File(fname1);

dfs_ref = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.smooth.' + hemi + '.dfs'))

segl = 225
sub1seg1 = normdata(fmri_sub1['func_'+hemi][:, :segl])
sub1seg2 = normdata(fmri_sub1['func_'+hemi][:, segl:2*segl])
sub2seg1 = normdata(fmri_sub2['func_'+hemi][:, :segl])
sub2seg2 = normdata(fmri_sub2['func_'+hemi][:, segl:2*segl])


rho_before = sp.sum(sub1seg1*sub1seg2, axis=1)/sub1seg1.shape[1]

dfs_ref = patch_color_attrib(dfs_ref, rho_before, clim=[0, .7])
view_patch_vtk(dfs_ref, azimuth=90, elevation=180, roll=90,
               outfile='before_seg1to2_1.png')
view_patch_vtk(dfs_ref, azimuth=-90, elevation=180, roll=-90,
               outfile='before_seg1to2_2.png')

_, Rot12 = rot_sub_data(ref=sub2seg2, sub=sub1seg1)

sub1seg1rot = sp.dot(sub1seg1, Rot12.T)

rho_after = sp.sum(sub1seg1rot*sub1seg2, axis=1)/sub2seg2.shape[1]
dfs_ref = patch_color_attrib(dfs_ref, rho_after, clim=[0, .7])
view_patch_vtk(dfs_ref, azimuth=90, elevation=180, roll=90,
               outfile='after_seg1to2_1.png')
view_patch_vtk(dfs_ref, azimuth=-90, elevation=180, roll=-90,
               outfile='after_seg1to2_2.png')

rho_direct22 = sp.sum(sub1seg2*sub2seg2, axis=1)/sub2seg2.shape[1]
dfs_ref = patch_color_attrib(dfs_ref, rho_direct22, clim=[0, .7])
sp.savez('movie_corr.npz', rho_direct22=rho_direct22)
view_patch_vtk(dfs_ref, azimuth=90, elevation=180, roll=90,
               outfile='direct_seg2to2_1.png')
view_patch_vtk(dfs_ref, azimuth=-90, elevation=180, roll=-90,
               outfile='direct_seg2to2_2.png')

rho_direct22_21 = rho_direct22 - rho_after
dfs_ref = patch_color_attrib(dfs_ref, rho_direct22_21, clim=[-.1, .1])
view_patch_vtk(dfs_ref, azimuth=90, elevation=180, roll=90,
               outfile='direct_seg2to2_1_sub.png')
view_patch_vtk(dfs_ref, azimuth=-90, elevation=180, roll=-90,
               outfile='direct_seg2to2_2_sub.png')
