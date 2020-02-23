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

p_dir = '/data_disk/HCP_data/data'
p_dir_ref = '/data_disk/HCP_data'
lst = os.listdir(p_dir)
r_factor = 3
ref_dir = os.path.join(p_dir_ref, 'reference')
nClusters = 2

ref = '196750'
sub = '196750'
print(sub, ref)
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = scipy.io.loadmat(fname1)  # h5py.File(fname1);
dfs_left = readdfs(
    os.path.join(p_dir_ref, 'reference', ref + '.aparc.a2009s.\
32k_fs.reduce3.left.dfs'))
dfs_left_sm = readdfs(
    os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.very_smooth.left.dfs'))

datasub = scipy.io.loadmat(
    os.path.join(p_dir, sub, sub + '.rfMRI_REST1_RL.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
dataref = scipy.io.loadmat(
    os.path.join(p_dir, ref, ref + '.rfMRI_REST2_RL.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))

LR_flag = msk['LR_flag']
LR_flag = np.squeeze(LR_flag) > 0
data = dataref['ftdata_NLM']
sub1, _, _ = normalizeData(data[LR_flag, :].T)

data = datasub['ftdata_NLM'][LR_flag, :]
#ind = sp.std(data, axis=1) > 0 #1e-116
#perm1 = sp.random.permutation(len(ind))
#data[ind, :] = data[ind[perm1], :]
sub2, _, _ = normalizeData(data.T)
sub2[500:, :] = 0
rho = np.sum(sub1 * sub2, axis=0)
dfs_left_sm = patch_color_attrib(dfs_left_sm, rho, clim=[-1, 1])
#dfs_left_sm.vColor[sp.absolute(rho) < 1e-116, :] = 0.5
view_patch_vtk(dfs_left_sm,
               azimuth=90,
               elevation=180,
               roll=90,
               outfile='sub1to2_view1.png',
               show=1)
view_patch_vtk(dfs_left_sm,
               azimuth=-90,
               elevation=-180,
               roll=-90,
               outfile='sub1to2_view2.png',
               show=1)

sub_rot, R1 = brainSync(X=sub1, Y=sub2)
rho = sp.sum(sub1 * sub_rot, axis=0)
dfs_left.attributes = rho
dfs_left_sm = patch_color_attrib(dfs_left_sm, rho, clim=[-1, 1])
dfs_left_sm.vColor[sp.absolute(rho) < 1e-116, :] = 0.5
print(sp.mean(rho[sp.absolute(rho) > 1e-116]))
view_patch_vtk(dfs_left_sm,
               azimuth=90,
               elevation=180,
               roll=90,
               outfile='sub1to2_view1_rot_r2.png',
               show=1)
view_patch_vtk(dfs_left_sm,
               azimuth=-90,
               elevation=-180,
               roll=-90,
               outfile='sub1to2_view2_rot_r2.png',
               show=1)
