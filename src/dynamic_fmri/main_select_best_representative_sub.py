# ||AUM||
import scipy.io as spio
import scipy as sp
import numpy as np
from fmri_methods_sipi import rot_sub_data
from surfproc import view_patch_vtk, patch_color_attrib
from dfsio import readdfs
import os
from brainsync import normalizeData, brainSync
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import h5py
import csv
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




#%%
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
#sub_data_orig = sub_data.copy()
from joblib import Parallel, delayed
import multiprocessing
    
# what are your inputs, and what operation do you want to 
# perform on each input. For example...
ind2 = np.arange(nSub)
num_cores = multiprocessing.cpu_count()    
 
for ind1 in range(nSub):  
#    results = Parallel(n_jobs=2)(delayed(sp.linalg.norm)(sub_data[:,:,ind1]-sub_data[:,:,i]) for i in ind2)
#    dist_all_orig[ind1, :] = results
#%%    
    for ind2 in range(nSub):
        dist_all_orig[ind1, ind2] = sp.linalg.norm(sub_data[:, :, ind1] -
                                                   sub_data[:, :, ind2])
        sub_data_rot, _ = brainSync(X=sub_data[:, :, ind1],
                                    Y=sub_data[:, :, ind2])
        dist_all_rot[ind1, ind2] = sp.linalg.norm(sub_data[:, :, ind1] -
                                                  sub_data_rot)
        print(ind1, ind2)


sp.savez('ADHD_pairwise_dist.npz', dist_all_rot=dist_all_rot,
         dist_all_orig=dist_all_orig, lst=lst)
######
#%%
a = sp.load('ADHD_pairwise_dist_all_sub_by_sub2.npz')
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
diff = 0
for ind in range(nsub):
    Y2, _ = brainSync(X=sub_data[:, :, q], Y=sub_data[:, :, ind])
    diff += (Y2 - sub_data[:, :, q]) ** 2

spio.savemat('diff_individual_atlas.mat', {'diff': diff})

#%% Create Average atlas by synchronizing everyones data to one subject
atlas = 0
for ind in range(nsub):
    Y2, _ = brainSync(X=sub_data[:, :, q], Y=sub_data[:, :, ind])
    atlas += Y2
atlas /= nsub

diff = 0
for ind in range(nsub):
    Y2, _ = brainSync(X=atlas, Y=sub_data[:, :, ind])
    diff += (Y2 - atlas) ** 2
    print(ind,)

spio.savemat('diff_average_atlas.mat', {'diff': diff})


#%% Compute difference for the virtual subject
diff = 0
r = h5py.File('/big_disk/ajoshi/fmri_Atlas_Result/result_raw.mat')
atlas = sp.array(r['X2']).T
for ind in range(nsub):
    Y2, _ = brainSync(X=atlas, Y=sub_data[:, :, ind])
    diff += (Y2 - atlas) ** 2
spio.savemat('diff_avg_atlas.mat', {'diff': diff})
