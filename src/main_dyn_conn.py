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
from surfproc import view_patch_vtk, patch_color_attrib

BFPPATH = '/big_disk/ajoshi/coding_ground/bfp'
BrainSuitePath = '/home/ajoshi/BrainSuite17a/svreg'

NCMP = 21

surfObj = readdfs(join(BFPPATH, 'supp_data', 'bci32kleft.dfs'))
numVert = len(surfObj.vertices)

sub1n = '/deneb_disk/HCP/196750/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii'
sub1n_tsk = '/deneb_disk/HCP/196750/MNINonLinear/Results/tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_Atlas.dtseries.nii'

#sub1 = nilearn.image.load_img(sub1n)
sub1 = nib.cifti2.cifti2.load(sub1n)
X = sub1.get_data().T
X, _, _ = normalizeData(X.T)


sub1tsk = nib.cifti2.cifti2.load(sub1n_tsk)
Xtsk = sub1tsk.get_data().T
Xtsk, _, _ = normalizeData(Xtsk.T)

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

a = nib.cifti2.Cifti2Image(Xnew, sub1.header, file_map=sub1.file_map)
a.to_filename('outfile_task.nii')

a = nib.cifti2.Cifti2Image(Xtsk-Xnew, sub1.header, file_map=sub1.file_map)
a.to_filename('outfile_diff.nii')


#%%
fname1 = 'left_motor1.png'
fname2 = 'left_motor2.png'

p_dir_ref = '/big_disk/ajoshi/HCP_data/'
ref = '196750'  # '100307'

Xtsk, _, _ = normalizeData(Xtsk)
Xnew, _, _ = normalizeData(Xnew)


lsurf = readdfs(os.path.join(p_dir_ref, 'reference', ref + '.aparc\
.a2009s.32k_fs.very_smooth.left.dfs'))
lind = np.where(lsurf.labels > 0)[0]
lsurf.attributes = np.zeros((lsurf.vertices.shape[0]))
lsurf.attributes[lind] = np.linalg.norm(Xtsk[50:150, :len(lind)]-Xnew[50:150, :len(lind)], axis=0)
lsurf = patch_color_attrib(lsurf, clim=[.45, 1.0])
view_patch_vtk(lsurf, azimuth=90, elevation=180, roll=90,
               outfile=fname1, show=0)
view_patch_vtk(lsurf, azimuth=-90, elevation=180, roll=-90,
               outfile=fname2, show=0)


