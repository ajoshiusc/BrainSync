# ||AUM||
import scipy.io
import scipy as sp
import numpy as np
from fmri_methods_sipi import rot_sub_data
from surfproc import view_patch_vtk, patch_color_attrib
from dfsio import readdfs
import os

p_dir = '/big_disk/ajoshi/HCP_data/data'
p_dir_ref = '/big_disk/ajoshi/HCP_data'
lst = os.listdir(p_dir)

r_factor = 3
ref_dir = os.path.join(p_dir_ref, 'reference')
nClusters = 30

ref = '100307'
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = scipy.io.loadmat(fname1)  # h5py.File(fname1);
dfs_left = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.left.dfs'))
dfs_left_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.very_smooth.left.dfs'))
count1 = 0
rho_rho = []
rho_all = []
#lst=lst[:1]
labs_all = sp.zeros((len(dfs_left.labels), len(lst)))
lst = [lst[0]]
sub = lst[0]
data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_LR.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
LR_flag = msk['LR_flag']
LR_flag = np.squeeze(LR_flag) != 0
data = data['ftdata_NLM']
temp = data[LR_flag, :]
temp[1, :] = sp.randn(1, temp.shape[1]) # temp[1000, :]
m = np.mean(temp, 1)
temp = temp - m[:, None]
s = np.std(temp, 1)+1e-16
temp = temp/s[:, None]
d1 = temp

#data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_RL.\
#reduce3.ftdata.NLM_11N_hvar_25.mat'))
#LR_flag = msk['LR_flag']
#LR_flag = np.squeeze(LR_flag) != 0
#data = data['ftdata_NLM']
#temp = data[LR_flag, :]
#temp[1, :] = sp.randn(1, temp.shape[1]) # temp[1000, :]
#m = np.mean(temp, 1)
#temp = temp - m[:, None]
#s = np.std(temp, 1)+1e-16
temp = temp/s[:, None]
perm1 = np.random.permutation(temp.shape[1])
d2 = temp[:,perm1]

if count1 == 0:
    sub_data1 = sp.zeros((d1.shape[0], d1.shape[1], len(lst)))
    sub_data2 = sp.zeros((d2.shape[0], d2.shape[1], len(lst)))

sub_data1[:, :, count1] = d1
sub_data2[:, :, count1] = d2

count1 += 1
print count1,

nSub = sub_data1.shape[2]
dist_all_orig = sp.zeros(len(dfs_left_sm.vertices))
dist_all_rot = dist_all_orig.copy()
sub_data_orig1 = sub_data1.copy()
sub_data_orig2 = sub_data2.copy()

for ind in range(nSub):
    dist_all_orig += sp.mean((sub_data_orig1[:, :, ind]-sub_data_orig2
                             [:, :, ind])**2.0, axis=(1))
    sub_data2[:, :, ind], R = rot_sub_data(ref=sub_data1[:, :, ind],
                                        sub=sub_data2[:, :, ind])
    dist_all_rot += sp.mean((sub_data1[:, :, ind]-sub_data2[:, :, ind])**2.0,
                            axis=(1))
    print ind,

dist_all_rot = dist_all_rot/(nSub)
dist_all_orig = dist_all_orig/(nSub)

var_all = sp.zeros((sub_data1.shape[0], sub_data2.shape[1]))

avg_sub_data = sp.mean(sub_data1, axis=2)

#dfs_left_sm = patch_color_attrib(dfs_left_sm, dist_all_orig, clim=[0, 1])
#view_patch_vtk(dfs_left_sm, azimuth=-90, elevation=-180,
#               roll=-90, outfile='dist_sess_orig_view1_1sub_left_permuted.png', show=0)
#view_patch_vtk(dfs_left_sm, azimuth=90, elevation=180, roll=90,
#               outfile='dist_sess_orig_view2_1sub_left_permuted.png', show=0)
#
#dfs_left_sm = patch_color_attrib(dfs_left_sm, dist_all_rot, clim=[0, 1])
#view_patch_vtk(dfs_left_sm, azimuth=-90, elevation=-180, roll=-90,
#               outfile='dist_sess_rot_view1_1sub_left_permuted.png', show=0)
#view_patch_vtk(dfs_left_sm, azimuth=90, elevation=180, roll=90,
#               outfile='dist_sess_rot_view2_1sub_left_permuted.png', show=0)

ind=sp.where(R>0.0001)
diff1=ind[0]-perm1[ind[1]]
print 'diff between orig and recovered permutation is :',sp.linalg.norm(diff1)
