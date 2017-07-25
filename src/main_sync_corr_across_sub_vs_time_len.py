# ||AUM||
import scipy.io
import scipy as sp
import numpy as np
from fmri_methods_sipi import rot_sub_data
from surfproc import view_patch_vtk, patch_color_attrib
from dfsio import readdfs
import os
import matplotlib.pyplot as plt
from brainsync import normalizeData, brainSync


p_dir = '/big_disk/ajoshi/HCP_data/data'
p_dir_ref = '/big_disk/ajoshi/HCP_data'
lst = os.listdir(p_dir)

r_factor = 3
ref_dir = os.path.join(p_dir_ref, 'reference')
nClusters = 30

ref = '196750'
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = scipy.io.loadmat(fname1)  # h5py.File(fname1);
dfs_left = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.left.dfs'))
dfs_left_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.very_smooth.left.dfs'))

# view_patch_vtk(dfs_left_sm)
rho_rho = []
rho_all = []

subind = 0
IntV = range(10, 1210, 10)
rho = sp.zeros((len(lst)-1, len(IntV)))
rho_orig = sp.zeros((len(lst)-1, len(IntV)))

for subno in range(len(lst)-1):
    sub = lst[subno]
    data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_LR.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
    LR_flag = msk['LR_flag']
    LR_flag = np.squeeze(LR_flag) != 0
    data = data['ftdata_NLM']
    temp = data[LR_flag, :]

    d1 = temp.T
    sub = lst[subno+1]
    data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_LR.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
    LR_flag = msk['LR_flag']
    LR_flag = np.squeeze(LR_flag) != 0
    data = data['ftdata_NLM']
    temp = data[LR_flag, :]

    d2 = temp.T

    ind = 0
    for len1 in IntV:
        sub_data1, _, _ = normalizeData(d1[:len1, :])
        sub_data2, _, _ = normalizeData(d2[:len1, :])
        s = sp.std(sub_data2, axis=0)
        sub_data1 = sub_data1[:, s > 1e-2]
        sub_data2 = sub_data2[:, s > 1e-2]
        sub_data2_sync, Rot = brainSync(X=sub_data1, Y=sub_data2)
        rho_orig[subno, ind] = sp.mean(sp.sum(sub_data2*sub_data1, axis=0))
        rho[subno, ind] = sp.mean(sp.sum(sub_data2_sync*sub_data1, axis=0))
        print(subno, len1, rho_orig[subno, ind], rho[subno, ind])
        ind += 1

    sp.savez('avg_corr_sub.npz', rho=rho, rho_orig=rho_orig)

plt.plot(IntV, rho)
plt.ylim(ymax=1, ymin=0.5)
plt.savefig('rho_sync_vs_len_diff_sub2.pdf')
plt.show()
