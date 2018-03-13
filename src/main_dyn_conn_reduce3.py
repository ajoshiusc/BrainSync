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

NCMP = 51

surfObj = readdfs(join(BFPPATH, 'supp_data', 'bci32kright.dfs'))
numVert = len(surfObj.vertices)

#sub1n='/big_disk/ajoshi/HCP100/HCP100/135932/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii';
sub1n='/big_disk/ajoshi/with_andrew/100307/100307.rfMRI_REST1_LR.reduce3.ftdata.hvar_0.mat'
sub1n_tsk='/big_disk/ajoshi/with_andrew/100307/100307.tfMRI_MOTOR_LR.reduce3.ftdata.hvar_0.mat'
#sub1n = '/deneb_disk/HCP/196750/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii'
#sub1n_tsk = '/deneb_disk/HCP/196750/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas.dtseries.nii'
X = spio.loadmat(sub1n)
X = X['ftdata']
X, _, _ = normalizeData(X.T)

#sub1 = nilearn.image.load_img(sub1n)
#sub1 = nib.load(sub1n)
#X = sub1.get_data().T

Xtsk = spio.loadmat(sub1n_tsk)
Xtsk = Xtsk['ftdata']
Xtsk, _, _ = normalizeData(Xtsk.T)

#sub1tsk = nib.cifti2.cifti2.load(sub1n_tsk)


# %% Explained variance
_, s, _ = np.linalg.svd(np.dot(X, X.T))
plt.figure()
plt.plot(s[:50])
plt.title('sigma plot')

# %% Do the PCA
p = PCA(n_components=NCMP)
D = p.fit_transform(X.T)

print("Explained Variance Fraction = %f" % p.explained_variance_ratio_.sum())

# D is the exeplar data
D, _, _ = normalizeData(D.T)

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

Xtsk, _, _ = normalizeData(Xtsk)
Xnew, _, _ = normalizeData(Xnew)

lsurf = surfObj 
ls = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.reduce3.very_smooth.right.dfs'))
lsurf=ls;
#lind = np.where(ls.labels > -10)[0]
lsurf.attributes = np.zeros((lsurf.vertices.shape[0]))
#lsurf.attributes = X[150,:lsurf.vertices.shape[0]] # 
lsurf.attributes = np.sum((Xtsk-Xnew)**2, axis=0)
lsurf.attributes = lsurf.attributes[lsurf.vertices.shape[0]:]
lsurf.attributes = smooth_surf_function(lsurf, lsurf.attributes)#, a1=1.1, a2=1.1)
lsurf = patch_color_attrib(lsurf, clim=[1, 2])
view_patch_vtk(lsurf, azimuth=90, elevation=180, roll=90,
               outfile=fname1, show=0)
view_patch_vtk(lsurf, azimuth=-90, elevation=180, roll=-90,
               outfile=fname2, show=0)


