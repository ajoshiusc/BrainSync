# ||AUM||
import scipy.io as spio
import scipy as sp
import numpy as np
from fmri_methods_sipi import rot_sub_data
from surfproc import view_patch_vtk, patch_color_attrib, smooth_surf_function, smooth_patch
from dfsio import readdfs
import os
from brainsync import normalizeData, brainSync
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import h5py
import csv
BFPPATH = '/big_disk/ajoshi/coding_ground/bfp'
BrainSuitePath = '/home/ajoshi/BrainSuite17a/svreg'

#%%

p_dir = '/deneb_disk/ADHD_Peking_bfp'
lst = os.listdir(p_dir)
count1 = 0
nsub = 0

#%% Read CSV File
normSub = [];adhdCombinedSub=[];adhdHyperactiveSub=[];adhdInattentive=[];
with open('/deneb_disk/ADHD_Peking_bfp/Peking_1_phenotypic.csv', newline='') as csvfile:    
    creader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for row in creader:
        dx=row['DX']
        sub=row['ScanDir ID']
        qc=row['QC_Rest_1']
        fname = os.path.join(p_dir, sub + '_rest_bold.32k.GOrd.mat')        

        if not os.path.isfile(fname) or int(qc) != 1:
            continue

        if int(dx) == 0:
            normSub.append(sub)

        if int(dx) == 1:
            adhdCombinedSub.append(sub)
            
        if int(dx) == 2:
            adhdHyperactiveSub.append(sub)
            
        if int(dx) == 3:
            adhdInattentive.append(sub)

        print(sub, dx, qc)




#%% Read Normal Subjects
count1 = 0
for sub in normSub:
    fname = os.path.join(p_dir, sub + '_rest_bold.32k.GOrd.mat')
    df = spio.loadmat(fname)
    data = df['dtseries'].T
    d, _, _ = normalizeData(data)

    if count1 == 0:
        sub_data = sp.zeros((d.shape[0], d.shape[1], len(normSub)))

    sub_data[:, :, count1] = d
    count1 += 1
    print(count1, )

#%%
nSub = sub_data.shape[2]
print(nSub)
dist_all_orig = sp.zeros([nSub, nSub])
dist_all_rot = dist_all_orig.copy()
 
for ind1 in range(nSub):  
    for ind2 in range(nSub):
        dist_all_orig[ind1, ind2] = sp.linalg.norm(sub_data[:, :, ind1] -
                                                   sub_data[:, :, ind2])
        sub_data_rot, _ = brainSync(X=sub_data[:, :, ind1],
                                    Y=sub_data[:, :, ind2])
        dist_all_rot[ind1, ind2] = sp.linalg.norm(sub_data[:, :, ind1] -
                                                  sub_data_rot)
        print(ind1, ind2)


sp.savez('ADHD_pairwise_dist.npz', dist_all_rot=dist_all_rot,
         dist_all_orig=dist_all_orig, normSub=normSub)
######
#%%
a = sp.load('ADHD_pairwise_dist.npz')
lst = a['lst']
q = sp.argmin(a['dist_all_rot'][:-1, :-1].sum(1))
print('The representative subject is: %s ' % lst[q])
m = MDS(n_components=2, dissimilarity='precomputed')
e = m.fit_transform(a['dist_all_rot'])
print(e)
fig, ax = plt.subplots()
ax.scatter(e[:, 0], e[:, 1])
for i in range(e.shape[0]):
    ax.annotate(lst[i], (e[i, 0], e[i, 1]))

#%% Compute difference
diff = 0; q=3
for ind in range(15):
    Y2, _ = brainSync(X=sub_data[:, :, q], Y=sub_data[:, :, ind])
    diff += (Y2 - sub_data[:, :, q]) ** 2
    print(ind, end=',')

spio.savemat('ADHD_norm_diff2sub1.mat', {'diff': diff})

#%% Create Average atlas by synchronizing everyones data to one subject
atlas = 0
nSub = len(normSub)
for ind in range(nSub-15):
    Y2, _ = brainSync(X=sub_data[:, :, q], Y=sub_data[:, :, ind])
    atlas += Y2
atlas /= (nSub-15)
spio.savemat('ADHD_avg_atlas.mat', {'atlas':atlas})
#%% Atlas to normal subjects diff
diff = 0
for ind in range(nSub-15,nSub):
    Y2, _ = brainSync(X=atlas, Y=sub_data[:, :, ind])
    diff += (Y2 - atlas) ** 2
    print(ind,)

spio.savemat('ADHD_diff_avg_atlas.mat', {'diff': diff})

#%% Read ADHD Inattentive
#del sub_data
count1 = 0
for sub in adhdInattentive:
    fname = os.path.join(p_dir, sub + '_rest_bold.32k.GOrd.mat')
    df = spio.loadmat(fname)
    data = df['dtseries'].T
    d, _, _ = normalizeData(data)
    if count1 == 0:
        sub_data = sp.zeros((d.shape[0], d.shape[1], len(adhdInattentive)))
        
    sub_data[:, :, count1] = d
    count1 += 1
    print(count1, )
#%% Atlas to normal subjects diff
diffAdhdInatt = 0

for ind in range(sub_data.shape[2]):
    Y2, _ = brainSync(X=atlas, Y=sub_data[:, :, ind])
    diffAdhdInatt += (Y2 - atlas) ** 2
    print(ind,)

spio.savemat('ADHD_diff_adhd_inattentive.mat', {'diffAdhdInatt': diffAdhdInatt})

#%% Read surfaces for visualization

lsurf = readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kleft.dfs')
rsurf = readdfs('/home/ajoshi/coding_ground/bfp/supp_data/bci32kright.dfs')
a=spio.loadmat('/home/ajoshi/coding_ground/bfp/supp_data/USCBrain_grayord_labels.mat')
labs=a['labels']
lsurf.attributes = np.zeros((lsurf.vertices.shape[0]))
rsurf.attributes = np.zeros((rsurf.vertices.shape[0]))
lsurf=smooth_patch(lsurf,iterations=3000)
rsurf=smooth_patch(rsurf,iterations=3000)
labs[sp.isnan(labs)]=0
diff=diff*(labs>0)
diffAdhdInatt=diffAdhdInatt*(labs>0)

nVert = lsurf.vertices.shape[0]

#%% Visualization of normal diff from the atlas
lsurf.attributes = np.sqrt(np.sum((diff), axis=0))
lsurf.attributes = lsurf.attributes[:nVert]/15
rsurf.attributes = np.sqrt(np.sum((diff), axis=0))
rsurf.attributes = rsurf.attributes[nVert:2*nVert]/15
lsurf = patch_color_attrib(lsurf, clim=[0.1,.3])
rsurf = patch_color_attrib(rsurf, clim=[0.1,.3])

view_patch_vtk(lsurf, azimuth=90, elevation=180, roll=90,
               outfile='l1normal.png', show=1)
view_patch_vtk(rsurf, azimuth=-90, elevation=180, roll=-90,
               outfile='r1normal.png', show=1)

#%% Visualization of ADHD diff from the atlas
lsurf.attributes = np.sqrt(np.sum((diffAdhdInatt), axis=0))
lsurf.attributes = lsurf.attributes[:nVert]/15
rsurf.attributes = np.sqrt(np.sum((diffAdhdInatt), axis=0))
rsurf.attributes = rsurf.attributes[nVert:2*nVert]/15
lsurf = patch_color_attrib(lsurf, clim=[0.1, .3])
rsurf = patch_color_attrib(rsurf, clim=[0.1, .3])

view_patch_vtk(lsurf, azimuth=90, elevation=180, roll=90,
               outfile='l1adhd.png', show=1)
view_patch_vtk(rsurf, azimuth=-90, elevation=180, roll=-90,
               outfile='r1adhd.png', show=1)

#%%
lsurf.attributes = np.sqrt(np.sum((diffAdhdInatt), axis=0))-np.sqrt(np.sum((diff), axis=0))
rsurf.attributes = np.sqrt(np.sum((diffAdhdInatt), axis=0))-np.sqrt(np.sum((diff), axis=0))
lsurf.attributes = lsurf.attributes[:nVert]/15
rsurf.attributes = rsurf.attributes[nVert:2*nVert]/15

lsurf.attributes = smooth_surf_function(lsurf,lsurf.attributes,3,3)
rsurf.attributes = smooth_surf_function(rsurf,rsurf.attributes,3,3)
lsurf = patch_color_attrib(lsurf, clim=[-0.01, 0.01])
rsurf = patch_color_attrib(rsurf, clim=[-0.01, 0.01])

view_patch_vtk(lsurf, azimuth=90, elevation=180, roll=90,
               outfile='l1adhd_normal_diff.png', show=1)
view_patch_vtk(rsurf, azimuth=-90, elevation=180, roll=-90,
               outfile='r1adhd_normal_diff.png', show=1)

#%%
pv = sp.zeros(diff.shape[1])
for vind in range(diff.shape[1]):
    _, pv[vind] = sp.stats.ranksums(diff[:, vind], diffAdhdInatt[:, vind])


lsurf.attributes = pv
rsurf.attributes = pv
lsurf.attributes = lsurf.attributes[:nVert]
rsurf.attributes = rsurf.attributes[nVert:2*nVert]

lsurf = patch_color_attrib(lsurf, clim=[0, .01])
rsurf = patch_color_attrib(rsurf, clim=[0, .01])

view_patch_vtk(lsurf, azimuth=90, elevation=180, roll=90,
               outfile='l1adhd_normal_pval.png', show=1)
view_patch_vtk(rsurf, azimuth=-90, elevation=180, roll=-90,
               outfile='r1adhd_normal_pval.png', show=1)

