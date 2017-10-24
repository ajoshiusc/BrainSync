# ||AUM||
import scipy.io
import scipy as sp
import numpy as np
from brainsync import normalizeData, brainSync
from surfproc import view_patch_vtk, patch_color_attrib
from dfsio import readdfs
import os

# This code synchronizes two sessions of each subject using left hemisphere \
# data and uses the estimated transform to sync the right hemisphere data

p_dir = '/big_disk/ajoshi/HCP_data/data'
p_dir_ref = '/big_disk/ajoshi/HCP_data'
lst = os.listdir(p_dir)

r_factor = 3
ref_dir = os.path.join(p_dir_ref, 'reference')

ref = '196750'
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = scipy.io.loadmat(fname1)  # h5py.File(fname1);
LR_flag = msk['LR_flag']
LR_flag = np.squeeze(LR_flag) == 0

dfs_right = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.right.dfs'))
dfs_right_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.very_smooth.right.dfs'))
count1 = 0

rho_rho = []
rho_all = []
labs_all = sp.zeros((len(dfs_right.labels), len(lst)))
avgCorr = 0
# Read all the data
for sub in lst:
    data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_LR.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
    data = data['ftdata_NLM']
    sub1L, _, _ = normalizeData(data[~LR_flag, :].T)
    sub1R, _, _ = normalizeData(data[LR_flag, :].T)

    data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST2_LR.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
    data = data['ftdata_NLM']
    sub2L, _, _ = normalizeData(data[~LR_flag, :].T)
    sub2R, _, _ = normalizeData(data[LR_flag, :].T)
    _, R = brainSync(X=sub1L, Y=sub2L)
    avgCorr += sp.mean(sub1R*sp.dot(R, sub2R), axis=1)
    count1 += 1
    print count1,

nSub = sub1L.shape[2]

avgCorr = sp.zeros(len(dfs_right_sm.vertices))

for ind in range(nSub):
    print ind,

avgCorr = avgCorr/(nSub)

dfs_right_sm = patch_color_attrib(dfs_right_sm, avgCorr, clim=[0, 1])
view_patch_vtk(dfs_right_sm, azimuth=-90, elevation=-180,
               roll=-90, outfile='corr_LR_right1.png', show=0)
view_patch_vtk(dfs_right_sm, azimuth=90, elevation=180,
               roll=90, outfile='corr_LR_right2.png', show=0)