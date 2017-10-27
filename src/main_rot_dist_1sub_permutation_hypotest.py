# ||AUM||
import scipy.io
import scipy as sp
import numpy as np
from surfproc import view_patch_vtk, patch_color_attrib
from dfsio import readdfs
import os

from brainsync import brainSync, normalizeData

p_dir = '/big_disk/ajoshi/HCP_data/data'
p_dir_ref = '/big_disk/ajoshi/HCP_data'
lst = os.listdir(p_dir)

ref_dir = os.path.join(p_dir_ref, 'reference')
r_factor = 3
ref = '196750'
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = scipy.io.loadmat(fname1)  # h5py.File(fname1);
dfs_left = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.left.dfs'))
dfs_left_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.very_smooth.right.dfs'))
count1 = 0
rho_rho = []
rho_all = []
labs_all = sp.zeros((len(dfs_left.labels), len(lst)))
lst = [lst[0]]
sub = lst[0]
data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_LR.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
LR_flag = msk['LR_flag']
LR_flag = np.squeeze(LR_flag) != 0
data = data['ftdata_NLM']
d1, _, _ = normalizeData(data[LR_flag, :].T)

data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST2_LR.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
data = data['ftdata_NLM']
temp, _, _ = normalizeData(data[LR_flag, :].T)

null_corr = sp.zeros((len(dfs_left_sm.vertices), 1000))

for iter1 in sp.arange(1000):
    perm1 = np.random.permutation(temp.shape[1])
    d2 = temp[:, perm1]

    d2, R = brainSync(X=d1, Y=d2)
    null_corr[:, iter1] = sp.sum(d1*d2, axis=0)
    print iter1,

d2, R = brainSync(X=d1, Y=temp)

scorr = sp.sum(d1*d2, axis=0)

c = scorr[:, None] < null_corr

pval = sp.mean(c, axis=1)

dfs_left_sm = patch_color_attrib(dfs_left_sm, 1-pval, clim=[0.95,1])
view_patch_vtk(dfs_left_sm, azimuth=-90, elevation=-180,
               roll=-90, outfile='pval_perm1.png', show=0)
view_patch_vtk(dfs_left_sm, azimuth=90, elevation=180, roll=90,
               outfile='pval_perm2.png', show=0)
