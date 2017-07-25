# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 05:39:43 2016

@author: ajoshi
"""
import scipy.io
import scipy as sp
import os
import numpy as np
from dfsio import readdfs
from surfproc import view_patch_vtk, patch_color_attrib
from fmri_methods_sipi import rot_sub_data
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from hrf import fast_fslgamma, double_gamma_hrf

p_dir = '/big_disk/ajoshi/HCP_data'
p_dir_ref='/big_disk/ajoshi/HCP_data/'
lst = os.listdir(p_dir)
r_factor = 3
ref_dir = os.path.join(p_dir_ref, 'reference')
nClusters = 3

ref = '196750'
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = scipy.io.loadmat(fname1)  # h5py.File(fname1);
dfs_left = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.left.dfs'))
dfs_left_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.very_smooth.left.dfs'))

ind_subsample = sp.arange(start=0, stop=dfs_left.labels.shape[0],step=100)
ind_rois_orig = sp.in1d(dfs_left.labels,[46,3,4,28,29,68,69,70])
ind_rois = sp.full(ind_rois_orig.shape[0], False ,dtype=bool)
ind_rois = ind_rois_orig.copy()
#ind_rois[ind_subsample] = True
hrf = double_gamma_hrf(sp.arange(50)*0.72) #fast_fslgamma
ind_rois = sp.nonzero(ind_rois)[0]

surf1 = dfs_left_sm

surf1.attributes = sp.zeros(surf1.vertices.shape[0])
surf1.attributes[ind_rois] = 1
surf1 = patch_color_attrib(surf1)
view_patch_vtk(surf1, show=1, azimuth=90, elevation=180, roll=90, outfile='motor_region.png')

# sub = '110411'
# p_dir = '/home/ajoshi/data/HCP_data'
lst = os.listdir('/big_disk/ajoshi/HCP5')
rho1 = 0; rho1rot = 0; rho2 = 0; rho2rot = 0;
# lst = [lst[0]]
diffbefore = 0
diffafter = 0

sub = lst[0]
vtscore = scipy.io.loadmat('/big_disk/ajoshi/with_andrew/100307/100307.reduce3.MOTOR_x100307_tfmri_motor_level2_t_avg_hp200_s4_tskdata.mat')

vrest = scipy.io.loadmat('/big_disk/ajoshi/with_andrew/100307/100307.tfMRI_MOTOR_RL.reduce3.ftdata.NLM_11N_hvar_5.mat')
#vrest = scipy.io.loadmat('/big_disk/ajoshi/with_andrew/100307/100307.tfMRI_LANGUAGE_RL.reduce3.ftdata.NLM_11N_hvar_0.45.mat')
#vrest = scipy.io.loadmat('//big_disk/ajoshi/with_andrew/100307/100307.rfMRI_REST2_LR.reduce3.ftdata.NLM_11N_hvar_5.mat')

LR_flag = msk['LR_flag']
LR_flag = np.squeeze(LR_flag) > 0
data = vrest['ftdata_NLM']
vrest = data[LR_flag]
vrest = vrest[ind_rois,]
m = np.mean(vrest, 1)
vrest = vrest - m[:, None]
s = np.std(vrest, 1)+1e-116
vmotor1 = vrest/s[:, None]

tval = sp.squeeze(vtscore['ftdata'][LR_flag])
tval[sp.isnan(tval)] = 0
surf1.attributes = sp.zeros(surf1.vertices.shape[0])
surf1.attributes[ind_rois] = tval[ind_rois]

ind_max_t = sp.argmax(surf1.attributes[ind_rois])
#ind_max_t = ind_rois[ind_max_t]
surf1 = patch_color_attrib(surf1)
view_patch_vtk(surf1, show=1, azimuth=90, elevation=180, roll=90, outfile='motor_tscore.png')

#vmotor1 = vmotor1[ind_rois,]
#vrest = nib.load('/big_disk/ajoshi/HCP5/' + sub + '/MNINonLinear/Resu\
#lts/rfMRI_REST2_LR/rfMRI_REST2_LR_Atlas_hp2000_clean.dtseries.nii')
vrest = scipy.io.loadmat('/big_disk/ajoshi/with_andrew/100307/100307.tfMRI_MOTOR_LR.reduce3.ftdata.NLM_11N_hvar_5.mat')
LR_flag = msk['LR_flag']
LR_flag = np.squeeze(LR_flag) > 0
data = vrest['ftdata_NLM']
#data = sp.squeeze(vrest.get_data()).T
vrest = data[LR_flag]
vrest = vrest[ind_rois,]
vrest = vrest[:, :vmotor1.shape[1]]    # make their length same
m = np.mean(vrest, 1)
vrest = vrest - m[:, None]
s = sp.std(vrest, axis=1) +1e-116
vmotor2 = vrest/s[:, None]

rho1 = sp.sum(vmotor2*vmotor1, axis=1)/vmotor2.shape[1]
diffbefore = vmotor2 - vmotor1

vmotor1orig=vmotor1.copy()
vmotor1, Rot = rot_sub_data(ref=vmotor2, sub=vmotor1) #, area_weight=sp.sqrt(surf_weight[ind_rois]))
rho1rot = sp.sum(vmotor2*vmotor1, axis=1)/vmotor2.shape[1]

diffafter = vmotor2 - vmotor1

#diffbefore = gaussian_filter(diffbefore,[0,5])
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
im=ax1.imshow(sp.absolute(diffbefore), aspect='auto', clim=(0, 2.0))
fig1.colorbar(im,ax=ax1)
fig1.savefig('dist_motor_before.png', dpi=300)

diffafter = gaussian_filter(diffafter, [0, 5])
fig2 = plt.figure()
ax1 = fig2.add_subplot(111)
im=ax1.imshow(sp.absolute(diffafter), aspect='auto', clim=(0, 2.0))
fig2.colorbar(im,ax=ax1)
fig2.savefig('dist_motor_after.png', dpi=300)

rho_full=sp.zeros((surf1.attributes.shape[0]))
rho_full[ind_rois] = rho1
dfs_left_sm.attributes = rho_full;
dfs_left_sm=patch_color_attrib(dfs_left_sm,clim=[0,1])
view_patch_vtk(dfs_left_sm, azimuth=90, elevation=180, roll=90, outfile='motor2motor_before_rot.png', show=1)

   #dfs_left_sm.attributes = sp.absolute(diffafter[:,t])
#    dfs_left_sm=patch_color_attrib(dfs_left_sm,clim=[0,1])
#    view_patch_vtk(dfs_left_sm, azimuth=90, elevation=180, roll=90, outfile='rest1motrho_full=sp.zeros((surf1.attributes.shape[0]))
rho_full[ind_rois] = rho1rot
dfs_left_sm.attributes = rho_full;
dfs_left_sm=patch_color_attrib(dfs_left_sm, clim=[0,1])
view_patch_vtk(dfs_left_sm, azimuth=90, elevation=180, roll=90, outfile='motor2motor_after_rot.png', show=1)
#view_patch_vtk(dfs_left_sm, azimuth=-90, elevation=180, roll=-90, outfile='hand_vs_hand_after_rot2.png', show=1)

aaa = vmotor1[ind_max_t, :]
aaaorig = vmotor1orig[ind_max_t, :]

bbb = vmotor2[ind_max_t, :]

import matplotlib.pyplot as plt

plt.plot(aaa,'r')
plt.plot(aaaorig,'k')
plt.plot(bbb,'b')
plt.show()

# Tongue
tblock1 = sp.zeros(vmotor1.shape[1])
tblock2 = sp.zeros(vmotor2.shape[1])
tblock1[57:74] = 1
tblock1[140:157] = 1

tblock2[78:95] = 1
tblock2[183:200] = 1

m = np.mean(tblock1)
tblock1 = tblock1 - m
s = sp.std(tblock1) +1e-116
tblock1 = tblock1/s

m = np.mean(tblock2)
tblock2 = tblock2 - m
s = sp.std(tblock2) +1e-116
tblock2 = tblock2/s

fig3 = plt.figure()
ax1 = fig3.add_subplot(111)

#blksynced = gaussian_filter(sp.dot(Rot,tblock2),[2])
blksynced = sp.convolve(hrf, sp.dot(Rot, tblock2), mode='same')
tblock2 = sp.convolve(hrf, tblock2, mode='same')
tblock1 = sp.convolve(hrf, tblock1, mode='same')

ax1.plot(tblock2,'r')
ax1.plot(tblock1,'k')

#blksynced = sp.dot(Rot,tblock2)
ax1.plot(blksynced*1.3,'b')
fig3.savefig('synced_blocks_tongue2.png',dpi=300)

# Right foot
tblock1 = sp.zeros(vmotor1.shape[1])
tblock2 = sp.zeros(vmotor2.shape[1])
tblock1[78:95] = 1
tblock1[246:263] = 1

tblock2[36:53] = 1
tblock2[204:221] = 1

m = np.mean(tblock1)
tblock1 = tblock1 - m
s = sp.std(tblock1) +1e-116
tblock1 = tblock1/s

m = np.mean(tblock2)
tblock2 = tblock2 - m
s = sp.std(tblock2) +1e-116
tblock2 = tblock2/s
#tblock1=gaussian_filter(tblock1,10)
#tblock2=gaussian_filter(tblock2,10)
blksynced = sp.convolve(hrf,sp.dot(Rot,tblock2),mode='same')
tblock2 = sp.convolve(hrf,tblock2,mode='same')
tblock1 = sp.convolve(hrf,tblock1,mode='same')

fig4 = plt.figure()
ax1 = fig4.add_subplot(111)

ax1.plot(tblock2,'r')
ax1.plot(tblock1,'k')
#blksynced = gaussian_filter(sp.dot(Rot,tblock2),[2])

ax1.plot(blksynced*1.3,'b')
fig4.savefig('synced_blocks_right_foot2.png',dpi=300)

#righthand
tblock1 = sp.zeros(vmotor1.shape[1])
tblock2 = sp.zeros(vmotor2.shape[1])
tblock1[15:32] = 1
tblock1[182:199] = 1

tblock2[121:138] = 1
tblock2[225:242] = 1

m = np.mean(tblock1)
tblock1 = tblock1 - m
s = sp.std(tblock1) +1e-116
tblock1 = tblock1/s

m = np.mean(tblock2)
tblock2 = tblock2 - m
s = sp.std(tblock2) + 1e-116
tblock2 = tblock2/s
#tblock1=gaussian_filter(tblock1,10)
#tblock2=gaussian_filter(tblock2,10)
fig5 = plt.figure()
ax1 = fig5.add_subplot(111)

blksynced = sp.convolve(hrf, sp.dot(Rot, tblock2), mode='same')
tblock2 = sp.convolve(hrf, tblock2, mode='same')
tblock1 = sp.convolve(hrf, tblock1, mode='same')

ax1.plot(tblock2, 'r')
ax1.plot(tblock1, 'k')
#blksynced = gaussian_filter(sp.dot(Rot,tblock2),[2])

plt.plot(blksynced*1.3, 'b')
fig5.savefig('synced_blocks_right_hand2.png', dpi=300)
#plt.show()


#lefthand
tblock1 = sp.zeros(vmotor1.shape[1])
tblock2 = sp.zeros(vmotor2.shape[1])
tblock1[99:116] = 1
tblock1[225:242] = 1

tblock2[15:32] = 1
tblock2[162:179] = 1

m = np.mean(tblock1)
tblock1 = tblock1 - m
s = sp.std(tblock1) + 1e-116
tblock1 = tblock1/s

m = np.mean(tblock2)
tblock2 = tblock2 - m
s = sp.std(tblock2) + 1e-116
tblock2 = tblock2/s
#tblock1=gaussian_filter(tblock1,10)
#tblock2=gaussian_filter(tblock2,10)
fig6 = plt.figure()
ax1 = fig6.add_subplot(111)

ax1.plot(tblock2, 'r')
ax1.plot(tblock1, 'k')
#blksynced = gaussian_filter(sp.dot(Rot,tblock2),[2])

ax1.plot(blksynced*1.3, 'b')
fig6.savefig('synced_blocks_left_hand2.png', dpi=300)
#plt.show()
