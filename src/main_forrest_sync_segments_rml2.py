# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 05:39:43 2016

@author: ajoshi
"""
from fmri_methods_sipi import normdata, rot_sub_data
from scipy.io import loadmat
import scipy as sp
from dfsio import readdfs, writedfs
import os.path
from surfproc import view_patch_vtk, patch_color_attrib, view_patch

p_dir_ref = '/big_disk/ajoshi/HCP_data/'
hemi = 'right'
ref = '100307'

fmri_sub11 = loadmat('/deneb_disk/studyforrest/sub-01-run1\
/fmri_tnlm_5_reduce3_v2.mat')  # h5py.File(fname1);

fmri_sub21 = loadmat('/deneb_disk/studyforrest/sub-02-run1\
/fmri_tnlm_5_reduce3_v2.mat')  # h5py.File(fname1);

fmri_sub12 = loadmat('/deneb_disk/studyforrest/sub-01-run2\
/fmri_tnlm_5_reduce3_v2.mat')  # h5py.File(fname1);

fmri_sub22 = loadmat('/deneb_disk/studyforrest/sub-02-run2\
/fmri_tnlm_5_reduce3_v2.mat')  # h5py.File(fname1);

dfs_ref = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.smooth.' + hemi + '.dfs'))


sub1seg1 = normdata(fmri_sub11['func_'+hemi][:, :420])
sub1seg2 = normdata(fmri_sub12['func_'+hemi][:, :420])
sub2seg1 = normdata(fmri_sub21['func_'+hemi][:, :420])
sub2seg2 = normdata(fmri_sub22['func_'+hemi][:, :420])


rho_before = sp.sum(sub1seg1*sub1seg2, axis=1)/sub1seg1.shape[1]

dfs_ref = patch_color_attrib(dfs_ref, rho_before, clim=[0, .7])
view_patch_vtk(dfs_ref, azimuth=90, elevation=180, roll=90,
               outfile='before2_seg1to2_1.png')
view_patch_vtk(dfs_ref, azimuth=-90, elevation=180, roll=-90,
               outfile='before2_seg1to2_2.png')

_, Rot12 = rot_sub_data(ref=sub2seg2, sub=sub1seg1)

sub1seg1rot = sp.dot(sub1seg1, Rot12.T)

rho_after = sp.sum(sub1seg1rot*sub1seg2, axis=1)/sub2seg2.shape[1]
dfs_ref = patch_color_attrib(dfs_ref, rho_after, clim=[0, .7])
writedfs('temp.dfs', dfs_ref)
view_patch_vtk(dfs_ref, azimuth=90, elevation=180, roll=90,
               outfile='after2_seg1to2_1.png')
view_patch_vtk(dfs_ref, azimuth=-90, elevation=180, roll=-90,
               outfile='after2_seg1to2_2.png')
rho_direct22 = sp.sum(sub1seg2*sub2seg2, axis=1)/sub2seg2.shape[1]
dfs_ref = patch_color_attrib(dfs_ref, rho_direct22, clim=[0, .7])

view_patch_vtk(dfs_ref, azimuth=90, elevation=180, roll=90,
               outfile='direct2_seg2to2_1.png')
view_patch_vtk(dfs_ref, azimuth=-90, elevation=180, roll=-90,
               outfile='direct2_seg2to2_2.png')
