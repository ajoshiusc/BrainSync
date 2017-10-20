# ||AUM||
import scipy.io
import scipy as sp
import numpy as np
#from fmri_methods_sipi import rot_sub_data
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

ref = '196750'#'100307'
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
#lst=lst[:1]
labs_all = sp.zeros((len(dfs_left.labels), len(lst)))
sub = lst[0]
data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_LR.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
LR_flag = msk['LR_flag']
LR_flag = np.squeeze(LR_flag) != 0
data = data['ftdata_NLM']
rest = data[LR_flag, :]
#temp[5000:6000, 500:700] += sp.rand(1000, 200) # temp[1000, :]
#m = np.mean(temp, 1)
#temp = temp - m[:, None]
#s = np.std(temp, 1)+1e-16
#temp = temp/s[:, None]
#d1 = temp

task = scipy.io.loadmat('/big_disk/ajoshi/with_andrew/100307/100307.tfMRI_MOTOR_LR.reduce3.ftdata.NLM_11N_hvar_5.mat')
task = task['ftdata_NLM']
temp = task[LR_flag, :] # rest[:, :task.shape[1]]

d1, _, _ = normalizeData(temp.T)

dist_all_rot = 0

for ind in range(40):
    sub = lst[ind]
    data = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_LR.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))
    LR_flag = msk['LR_flag']
    LR_flag = np.squeeze(LR_flag) != 0
    data = data['ftdata_NLM']
    temp = data[LR_flag, :]
#    m = np.mean(temp, 1)
#    temp = temp - m[:, None]
#    s = np.std(temp, 1)+1e-16
#    temp = temp/s[:, None]
    temp = temp[:, :d1.shape[0]]
    d2, _, _ = normalizeData(temp.T)

#    d2 = temp

    sub_data1 = d1.T
    sub_data2 = d2.T

#    ind1 = s>1e-10

    dist_all_orig = sp.zeros(len(dfs_left_sm.vertices))
#    dist_all_rot = dist_all_orig.copy()
    sub_data_orig1 = sub_data1.copy()
    sub_data_orig2 = sub_data2.copy()

    dist_all_orig = sub_data_orig1-sub_data_orig2
    sub_data2, _ = brainSync(X=sub_data1.T, Y=sub_data2.T)
    sub_data2 = sub_data2.T
    # rot_sub_data(ref=sub_data1, sub=sub_data2)
    dist_all_rot += sub_data1-sub_data2
    print ind,

plt.figure()
plt.set_cmap('jet')
plt.imshow(sp.absolute(dist_all_orig), aspect='auto', clim=(0, 0.1))
plt.colorbar()
plt.savefig('dist_before_task.pdf', dpi=300)
plt.show()
plt.figure()

plt.set_cmap('jet')
plt.imshow(sp.absolute(dist_all_rot/40), aspect='auto', clim=(0, 0.1))
plt.colorbar()
plt.savefig('dist_after_task.pdf', dpi=300)
plt.show()
