# AUM
# Shree Ganeshaya Namaha
from dfsio import readdfs
from os.path import join
#import nilearn.image
import numpy as np
from brainsync import normalizeData, brainSync
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import nibabel as nib
import os
from surfproc import view_patch_vtk, patch_color_attrib, smooth_surf_function
import scipy.io as spio
BFPPATH = '/big_disk/ajoshi/coding_ground/bfp'
BrainSuitePath = '/home/ajoshi/BrainSuite17a/svreg'

NCMP = 21

ref = '196750'
p_dir = '/big_disk/ajoshi/HCP_data/data'
p_dir_ref = '/big_disk/ajoshi/HCP_data'
ref_dir = os.path.join(p_dir_ref, 'reference')
lst = os.listdir(p_dir)

fn1 = ref + '.reduce3' + '.LR_mask.mat'
fname1 = os.path.join(ref_dir, fn1)
msk = spio.loadmat(fname1)  # h5py.File(fname1);

refdata = spio.loadmat(os.path.join(p_dir, ref, ref + '.rfMRI_REST1_LR.reduce3.ftdata.hvar_0.mat'))

dr = refdata['ftdata'].T
meanData = 0
NSUB = len(lst)

for ind in range(NSUB):
    sub = lst[ind]
    data = spio.loadmat(os.path.join(p_dir, sub, sub + '.rfMRI_REST1_LR.\
reduce3.ftdata.hvar_0.mat'))
    LR_flag = msk['LR_flag']
    LR_flag = np.squeeze(LR_flag) != 0
    data = data['ftdata']
 #   temp = data[LR_flag, :]
#    m = np.mean(temp, 1)
#    temp = temp - m[:, None]
#    s = np.std(temp, 1)+1e-16
#    temp = temp/s[:, None]
#    temp = temp[:, :d1.shape[0]]
    d2, _, _ = normalizeData(data.T)    
    d2, _ = brainSync(dr, d2)
    meanData = meanData + d2
    print(ind,end=',')

# %% Do the PCA
np.savez('mean_data.npz', meanData=meanData)
p = PCA(n_components=NCMP)
D = p.fit_transform(meanData.T).T


# %% Explained variance
_, s, _ = np.linalg.svd(np.dot(meanData, meanData.T))
plt.figure()
plt.plot(s[:50])
plt.title('sigma plot')

# %%
print("Explained Variance Fraction = %f" % p.explained_variance_ratio_.sum())

# D is the exeplar data
D, _, _ = normalizeData(D)

task = spio.loadmat('/big_disk/ajoshi/with_andrew/100307/100307.tfMRI_MOTOR_LR.reduce3.ftdata.hvar_0.mat')
task = task['ftdata']

Xtsk, _, _ = normalizeData(task.T)
nT = Xtsk.shape[0]

Xnew = np.zeros(Xtsk.shape)
Cind = np.int((NCMP-1)/2+1)
for i in range(Xtsk.shape[0]-NCMP):
    xin = Xtsk[i:i+NCMP, :]
    xin, _, nrm = normalizeData(xin)
    dd, _ = brainSync(xin, D)
    dd = dd*nrm
    Xnew[Cind+i, :] = dd[Cind, :]
    print("%d" % i, end=',')

Xnew, _, _ = normalizeData(Xnew)

#a = nib.cifti2.Cifti2Image(Xnew, sub1.header, file_map=sub1.file_map)
#a.to_filename('outfile_task.nii')
#
#a = nib.cifti2.Cifti2Image(Xtsk-Xnew, sub1.header, file_map=sub1.file_map)
#a.to_filename('outfile_diff.nii')

#loading cifti files has indices garbled
#%%
fname1 = 'right_motor1.png'
fname2 = 'right_motor2.png'

p_dir_ref = '/big_disk/ajoshi/HCP_data/'
ref = '196750'  # '100307'


#lsurf = surfObj 
ls = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc.\
a2009s.32k_fs.reduce3.very_smooth.right.dfs'))

lsurf=ls;
#lind = np.where(ls.labels > -10)[0]
lsurf.attributes = np.zeros((lsurf.vertices.shape[0]))
#lsurf.attributes = X[150,:lsurf.vertices.shape[0]] # 

nVert = lsurf.vertices.shape[0]
diffafter = Xtsk - Xnew

lsurf.attributes = np.sum((diffafter)**2, axis=0)
lsurf.attributes = lsurf.attributes[nVert:]
#lsurf.attributes = smooth_surf_function(lsurf, lsurf.attributes)#, a1=1.1, a2=1.1)
lsurf = patch_color_attrib(lsurf, clim=[1, 2])
view_patch_vtk(lsurf, azimuth=90, elevation=180, roll=90,
               outfile=fname1, show=0)
view_patch_vtk(lsurf, azimuth=-90, elevation=180, roll=-90,
               outfile=fname2, show=0)


#%%
for ind in np.arange(Xtsk.shape[0]):
    lsurf.attributes = (diffafter[ind, nVert:])**2 #diffafter[ind, nVert:]
    fname1 = 'rest2motor_right_%d_d.png' % ind
    fname2 = 'rest2motor_right_%d_m.png' % ind
#    lsurf.attributes = smooth_surf_function(lsurf, lsurf.attributes, a1=1.1, a2=1.1)
    lsurf = patch_color_attrib(lsurf, clim=[0, .02])
    view_patch_vtk(lsurf, azimuth=90, elevation=180, roll=90,
                   outfile=fname1, show=0)
    view_patch_vtk(lsurf, azimuth=-90, elevation=180, roll=-90,
                   outfile=fname2, show=0)
    print(ind,)
#


np.savez('motor_diff_data.npz', diffafter=diffafter)

spio.savemat('motor_diff_data.mat',{'diffafter':diffafter})


#    d2 = temp

 