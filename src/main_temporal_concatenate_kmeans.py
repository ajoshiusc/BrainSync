# ||AUM||
import scipy.io
import scipy as sp
import numpy as np
from fmri_methods_sipi import rot_sub_data
from dfsio import readdfs
import os
from sklearn.cluster import KMeans

p_dir = '/big_disk/ajoshi/HCP_data/data'
p_dir_ref = '/big_disk/ajoshi/HCP_data'
lst = os.listdir(p_dir)
r_factor = 3
ref_dir = os.path.join(p_dir_ref, 'reference')
nClusters = 100

ref = '196750'  # chosen by taking smallest distance to every other subject
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = scipy.io.loadmat(fname1)  # h5py.File(fname1);
dfs_left = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.a2009s.\
32k_fs.reduce3.left.dfs'))
dfs_left_sm = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.a2009\
s.32k_fs.reduce3.very_smooth.left.dfs'))
count1 = 0
rho_rho = []
rho_all = []

labs_all = sp.zeros((len(dfs_left.labels), len(lst)))

for sub in lst:
    data1 = scipy.io.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_RL.\
reduce3.ftdata.NLM_11N_hvar_25.mat'))

    LR_flag = msk['LR_flag']
    LR_flag = np.squeeze(LR_flag) > 0
    data1 = data1['ftdata_NLM']

    temp = data1[LR_flag, :]
    m = np.mean(temp, 1)
    temp = temp - m[:, None]
    s = np.std(temp, 1)+1e-16
    temp = temp/s[:, None]
    d1 = temp

    if count1 == 0:
        sub_data1 = sp.zeros((d1.shape[0], d1.shape[1], len(lst)))

    sub_data1[:, :, count1] = d1

    count1 += 1
    print count1,

nSub = sub_data1.shape[2]

cat_data = sp.reshape(sub_data1, (sub_data1.shape[0],
                                  sub_data1.shape[1]*nSub), 'F')

print cat_data.shape

del sub_data1, d1, temp, data1


SC = KMeans(n_clusters=nClusters, random_state=5324, verbose=1)
labs_cat = SC.fit_predict(cat_data)

sp.savez_compressed('labs_concat_data_100_clusters', labs_cat=labs_cat,
                    lst=lst, nClusters=nClusters)
