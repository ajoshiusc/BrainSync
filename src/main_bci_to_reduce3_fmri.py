# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 15:51:16 2016

@author: ajoshi
"""
import sys
#sys.path.append('/home/ajoshi/coding_ground/cortical_parcellation/src/')

import nibabel.freesurfer.io as fsio
from surfproc import view_patch_vtk, patch_color_attrib
from dfsio import writedfs, readdfs
from nibabel.gifti.giftiio import read as gread
import os
from scipy.spatial import cKDTree
import scipy.io
import nibabel as nib
import glob

p_dir_ref = '/home/ajoshi/HCP_data'
ref_dir = os.path.join(p_dir_ref, 'reference')
ref = '100307'


def interpolate_fmri(fromsurf=[], tosurf=[]):
    ''' interpolate labels from surface to to surface'''
    tree = cKDTree(fromsurf.vertices)
    d, inds = tree.query(tosurf.vertices, k=1, p=2)
    tosurf.fmri = fromsurf.fmri[inds, :]
    return tosurf


class s:
    pass


class bci:
    pass


lst = glob.glob('/deneb_disk/studyforrest/sub-10*run3')

for sub in lst:
#    if os.path.exists(sub + '/fmrit_reduce3.mat'):
#        continue
    print sub
    if not os.path.exists(sub + '/anat/BrainSuite/fmri_surf_dat_v2.mat'):
        continue

    ''' Right Hemisphere '''
    ''' BCI to FS processed BCI '''
    bci_BrainSuite = readdfs('/home/ajoshi/coding_ground/svreg/BCI-DNI\
_brain_atlas/BCI-DNI_brain.right.inner.cortex.dfs')

    data1 = scipy.io.loadmat(sub + '/anat/BrainSuite/fmri_surf_dat_v2.mat')
    bci_BrainSuite.fmri = data1['datar_atlas'].squeeze()
    bci_BrainSuite.vertices[:, 0] -= 96*0.8
    bci_BrainSuite.vertices[:, 1] -= 192*0.546875
    bci_BrainSuite.vertices[:, 2] -= 192*0.546875
    bci.vertices, bci.faces = fsio.read_geometry('/big_disk/ajoshi/data/BCI_\
DNI_Atlas/surf/rh.white')
    bci = interpolate_fmri(bci_BrainSuite, bci)

    ''' FS_BCI to FS BCI Sphere'''
    bci.vertices, bci.faces = fsio.read_geometry('/big_disk/ajoshi/data/BCI_\
DNI_Atlas/surf/rh.sphere.reg')
#    bci.attributes = bci.vertices[:,0]
#    bci = patch_color_attrib(bci)
#    view_patch_vtk(bci)
    ''' FS BCI Sphere to ref FS Sphere'''
    g_surf = nib.load('/big_disk/ajoshi/HCP_data/reference/100307/MNINon\
Linear/Native/100307.R.sphere.reg.native.surf.gii')
    s.vertices = g_surf.darrays[0].data
    s.faces = g_surf.darrays[1].data
    s = interpolate_fmri(bci, s)

    ''' ref BCI Sphere to FS very inflated '''
    g_surf = nib.load('/big_disk/ajoshi/HCP_data/reference/100307/MNINon\
Linear/Native/100307.R.very_inflated.native.surf.gii')
    bci.vertices = g_surf.darrays[0].data
    bci.faces = g_surf.darrays[1].data
    bci.fmri = s.fmri
    ''' FS very inflated to reduce3 '''
    dfs = readdfs('/big_disk/ajoshi/HCP_data/reference/100307.aparc.a\
2009s.32k_fs.reduce3.very_smooth.right.dfs')
    dfs = interpolate_fmri(bci, dfs)
    fmri = dfs.fmri
    a = {}
    a['fmri_right'] = fmri


    ''' Left Hemisphere '''
    ''' BCI to FS processed BCI '''
    bci_BrainSuite = readdfs('/home/ajoshi/coding_ground/svreg/BCI-DNI\
_brain_atlas/BCI-DNI_brain.left.inner.cortex.dfs')

#    data1 = scipy.io.loadmat(os.path.join(sub + '/anat/BrainSuite/fmri_surf_dat.mat'))
    bci_BrainSuite.fmri = data1['datal_atlas'].squeeze()
    bci_BrainSuite.vertices[:, 0] -= 96*0.8
    bci_BrainSuite.vertices[:, 1] -= 192*0.546875
    bci_BrainSuite.vertices[:, 2] -= 192*0.546875
    bci.vertices, bci.faces = fsio.read_geometry('/big_disk/ajoshi/data/BCI_\
DNI_Atlas/surf/lh.white')
    bci = interpolate_fmri(bci_BrainSuite, bci)

    ''' FS_BCI to FS BCI Sphere'''
    bci.vertices, bci.faces = fsio.read_geometry('/big_disk/ajoshi/data/BCI_\
DNI_Atlas/surf/lh.sphere.reg')

    ''' FS BCI Sphere to ref FS Sphere'''
    g_surf = nib.load('/big_disk/ajoshi/HCP_data/reference/100307/MNINon\
Linear/Native/100307.L.sphere.reg.native.surf.gii')
    s.vertices = g_surf.darrays[0].data
    s.faces = g_surf.darrays[1].data
    s = interpolate_fmri(bci, s)

    ''' ref BCI Sphere to FS very inflated '''
    g_surf = nib.load('/big_disk/ajoshi/HCP_data/reference/100307/MNINon\
Linear/Native/100307.L.very_inflated.native.surf.gii')
    bci.vertices = g_surf.darrays[0].data
    bci.faces = g_surf.darrays[1].data
    bci.fmri = s.fmri
    ''' FS very inflated to reduce3 '''
    dfs = readdfs('/big_disk/ajoshi/HCP_data/reference/100307.aparc.a\
2009s.32k_fs.reduce3.very_smooth.left.dfs')
    dfs = interpolate_fmri(bci, dfs)
    fmri = dfs.fmri
#    a = {}
    a['fmri_left'] = fmri
    scipy.io.savemat(sub + '/fmrit_reduce3_v2.mat', a)
    print(sub)
