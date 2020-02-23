# ||AUM||
from dfsio import readdfs
import scipy.io
import numpy as np
from mayavi import mlab
#import h5py
import os
#from scipy.stats import trim_mean
import scipy.sparse as sp
from sklearn.cluster import SpectralClustering, DBSCAN
import matplotlib.pylab as plt
from tqdm import tqdm
from surfproc import view_patch, view_patch_vtk, get_cmap, smooth_patch, patch_color_attrib, smooth_surf_function, patch_color_labels

p_dir = '/data_disk/HCP_data'
lst = os.listdir(p_dir + '/data')
r_factor = 3
ref_dir = os.path.join(p_dir, 'reference.old')

ref = '100307'
print(ref + '.reduce' + str(r_factor) + '.LR_mask.mat')
fn1 = ref + '.reduce' + str(r_factor) + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = scipy.io.loadmat(fname1)  # h5py.File(fname1);
dfs_left = readdfs(
    os.path.join(p_dir, 'reference.old',
                 ref + '.aparc.a2009s.32k_fs.reduce3.left.dfs'))
dfs_left_sm = readdfs(
    os.path.join(p_dir, 'reference.old',
                 ref + '.aparc.a2009s.32k_fs.reduce3.very_smooth.left.dfs'))
count1 = 0
rho_rho = []
rho_all = np.zeros([0, 0])
for sub in tqdm(lst[:5]):
    try:
        data = scipy.io.loadmat(
            os.path.join(
                p_dir, 'data', sub,
                sub + '.rfMRI_REST2_RL.reduce3.ftdata.NLM_11N_hvar_25.mat'))
    except:
        continue

    LR_flag = msk['LR_flag']
    LR_flag = np.squeeze(LR_flag) > 0
    data = data['ftdata_NLM']
    temp = data[LR_flag, :]
    m = np.mean(temp, 1)
    temp = temp - m[:, None]
    s = np.std(temp, 1) + 1e-16
    temp = temp / s[:, None]
    msk_small_region = (dfs_left.labels == 46) | (dfs_left.labels == 28) | (
        dfs_left.labels == 29)  # % motor
    d = temp[msk_small_region, :]
    if count1 == 0:
        data_all = d
    else:
        data_all = np.hstack([data_all, d])

    count1 += 1
#    if count1 > 0:
#        break

B = np.corrcoef(data_all)
B[~np.isfinite(B)] = 0
B = np.abs(B)
plt.figure()
plt.imshow(np.abs(data_all))
plt.show()
plt.figure()
plt.imshow(np.abs(B))
plt.show()

for nClusters in [2, 3, 4, 5, 6]:
    SC = SpectralClustering(n_clusters=nClusters, affinity='precomputed')
    #SC=SpectralClustering(n_clusters=nClusters,assign_labels='discretize')
    labs = SC.fit_predict(B)

    r = dfs_left_sm
    r.labels = r.labels * 0
    r.labels[msk_small_region] = labs + 1

    dfs_left_sm = patch_color_labels(dfs_left_sm)
    #dfs_left_sm.vColor[sp.absolute(rho) < 1e-116, :] = 0.5
    filename = 'c' + str(nClusters) + 'labels_1.png'
    view_patch_vtk(dfs_left_sm,
                   azimuth=90,
                   elevation=180,
                   roll=90,
                   outfile=filename,
                   show=1)

    filename = 'c' + str(nClusters) + 'labels_2.png'
    view_patch_vtk(dfs_left_sm,
                   azimuth=-90,
                   elevation=-180,
                   roll=-90,
                   outfile=filename,
                   show=1)
