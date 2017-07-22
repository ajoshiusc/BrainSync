# ||AUM||
import scipy.io
import scipy as sp
import numpy as np
from dfsio import readdfs, writedfs
from mayavi import mlab
from fmri_methods_sipi import rot_sub_data
#import h5py
import os
from surfproc import view_patch, view_patch_vtk, get_cmap, smooth_patch, patch_color_attrib, smooth_surf_function
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import DictionaryLearning
from scipy.stats import trim_mean
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.utils.linear_assignment_ import linear_assignment
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from brainsync import brainSync, normalizeData

p_dir = '/big_disk/ajoshi/HCP_data/data'
p_dir_ref='/big_disk/ajoshi/HCP_data'
lst = os.listdir(p_dir)
r_factor = 3
ref_dir = os.path.join(p_dir_ref, 'reference')
nClusters=2

ref = '100307'
sub = '100307' #lst[25]
print sub, ref
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = scipy.io.loadmat(fname1)  # h5py.File(fname1);
dfs_left = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.a2009s.32k_fs.reduce3.left.dfs'))
dfs_left_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.a2009s.32k_fs.reduce3.very_smooth.left.dfs'))
count1 = 0
roilist=[30, 72, 9, 47] #pc
#ref=lst[11]
datasub = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST2_RL.reduce3.ftdata.NLM_11N_hvar_25.mat'))
dataref = scipy.io.loadmat(os.path.join(p_dir, ref, ref + '.rfMRI_REST1_RL.reduce3.ftdata.NLM_11N_hvar_25.mat'))

LR_flag = msk['LR_flag']
LR_flag = np.squeeze(LR_flag) > 0
data = dataref['ftdata_NLM']
sub1,_,_ = normalizeData(data[LR_flag, :].T)
#sub1 = sub1.T

data = datasub['ftdata_NLM']
sub2,_,_ = normalizeData(data[LR_flag, :].T)
#sub2 = sub2.T

msk_small_region = np.in1d(dfs_left.labels,roilist)
#    msk_small_region = (dfs_left.labels == 30) | (dfs_left.labels == 72) | (dfs_left.labels == 9) |  (dfs_left.labels == 47)  # % motor
d = sub1[:,msk_small_region]

ref_mean_pc = sp.mean(d,axis=1)
ref_mean_pc=ref_mean_pc-sp.mean(ref_mean_pc)
ref_mean_pc=ref_mean_pc/(sp.std(ref_mean_pc))


d = sub2[:,msk_small_region]

sub_mean_pc = sp.mean(d,axis=1)
sub_mean_pc=sub_mean_pc-sp.mean(sub_mean_pc)
sub_mean_pc=sub_mean_pc/(sp.std(sub_mean_pc))


rho = np.dot(ref_mean_pc,sub1)/ref_mean_pc.shape[0]
rho[~np.isfinite(rho)] = 0

simil_mtx=sp.pi/2.0 + sp.arcsin(rho)
#    simil_mtx=0.3*sp.ones(rho.shape)
#SC = SpectralClustering(n_clusters=nClusters, affinity='precomputed')
#labs_all = SC.fit_predict(simil_mtx)+1

#rho = smooth_surf_function(dfs_left_sm, rho)
dfs_left_sm.attributes = rho
dfs_left_sm=patch_color_attrib(dfs_left_sm, rho, clim=[-1,1])
view_patch_vtk(dfs_left_sm, azimuth=90,elevation=180, roll=90, outfile='sub1to1_view1_pc.png', show=1)
view_patch_vtk(dfs_left_sm, azimuth=-90,elevation=-180, roll=-90, outfile='sub1to1_view2_pc.png', show=1)

rho = np.dot(sub_mean_pc,sub2)/sub_mean_pc.shape[0]
rho[~np.isfinite(rho)] = 0

dfs_left_sm.attributes = rho
dfs_left_sm=patch_color_attrib(dfs_left_sm, rho, clim=[-1,1])
view_patch_vtk(dfs_left_sm, azimuth=90,elevation=180, roll=90, outfile='sub2to2_view1_pc.png', show=1)
view_patch_vtk(dfs_left_sm, azimuth=-90,elevation=-180, roll=-90, outfile='sub2to2_view2_pc.png', show=1)



rho = np.dot(ref_mean_pc,sub2)/ref_mean_pc.shape[0]
rho[~np.isfinite(rho)] = 0
#rho = smooth_surf_function(dfs_left_sm, rho)
dfs_left.attributes = rho
dfs_left_sm=patch_color_attrib(dfs_left_sm, rho, clim=[-1,1])
view_patch_vtk(dfs_left_sm, azimuth=90,elevation=180, roll=90, outfile='sub1to2_view1_pc.png', show=1)
view_patch_vtk(dfs_left_sm, azimuth=-90,elevation=-180, roll=-90, outfile='sub1to2_view2_pc.png', show=1)

#sm=smooth_patch(dfs_left,iter=1000)
#view_patch(sm)
#view_patch(dfs_left,rho)

sub_rot, R1 = brainSync(X=sub1, Y=sub2)
#sub_rot2, R2, C2 = rot_sub_data(sub1.T, sub2.T)
#sub_rot=sub_rot.T
rho = sp.dot(ref_mean_pc,sub_rot)/ref_mean_pc.shape[0]
#rho=rho[0,1:]
rho[~np.isfinite(rho)] = 0
#rho = smooth_surf_function(dfs_left_sm, rho)
dfs_left.attributes = rho
dfs_left_sm=patch_color_attrib(dfs_left_sm, rho, clim=[-1,1])
view_patch_vtk(dfs_left_sm, azimuth=90,elevation=180, roll=90, outfile='sub1to2_view1_pc_rot.png', show=1)
view_patch_vtk(dfs_left_sm, azimuth=-90,elevation=-180, roll=-90, outfile='sub1to2_view2_pc_rot.png', show=1)


plt.figure(num=None, figsize=(10, 1), dpi=300, facecolor='w', edgecolor='k')
plt.plot(sub1[8991,100:200])
plt.plot(sub2[8991,100:200])
plt.plot(sub_rot[8991,100:200])
plt.savefig('Sync_time_series.pdf')
#plt.ylabel('time')
#plt.plot(sub1[8991,100:200], sub2[8991,100:200])
#plt.plot()
#plt.plot(sub1[8991,100:200],sub_rot[8991,100:200])
#plt.ylabel('time')
