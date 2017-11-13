# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 05:39:43 2016

@author: ajoshi
"""
import os
import glob
from dfsio import readdfs, writedfs
import scipy as sp
from brainsync import brainSync, normalizeData
from statsmodels.stats.multitest import multipletests
from surfproc import view_patch_vtk, patch_color_attrib


BFPPATH = '/home/ajoshi/coding_ground/bfp/supp_data'

# Read Reference parcellation
refLeft = readdfs(os.path.join(BFPPATH, 'bci32kleft.dfs'))
refRight = readdfs(os.path.join(BFPPATH, 'bci32kright.dfs'))
#
nullsubDir = '/deneb_disk/Beijing_Zhang_bfp/'
lst = glob.glob(nullsubDir+'*LB40.mat')
nsub = 50#len(lst)
vsub = sp.io.loadmat('/deneb_disk/from_Todd_Constable_Epilepsy_Processed\
/sn7602/func/sn7602_rest_bold.32k.GOrd_LB40.mat')
vsub, _, _ = normalizeData(vsub['dtseries'].T)

print("There are %d subjects" % nsub)
print("Reading the subject data")
# %%
for ind1 in range(nsub):
    vrest = sp.io.loadmat(os.path.join(nullsubDir, lst[ind1]))
    vrest = vrest['dtseries'].T
    if ind1 == 0:
        vrest_subs = sp.zeros([vrest.shape[1], vrest.shape[0], nsub])

    vrest = vrest[:vsub.shape[0], :]
    vrest_subs[:, :, ind1], _, _ = normalizeData(vrest)
#    print(ind1, end=' ')
    print ind1,

# %%
#  Build Null Distribution
nVert = vrest_subs.shape[1]
rho_null = sp.zeros([nsub, nVert])
for ind1 in range(nsub):
    t = 0
    for ind2 in range(nsub):
        if ind1 == ind2:
            continue
        vrest1 = vrest_subs[:, :, ind1]
        vrest2 = vrest_subs[:, :, ind2]
        vrest2, Rot = brainSync(X=vrest1, Y=vrest2)
        t += sp.sum(vrest1*vrest2, axis=0)

    rho_null[ind1, :] = t/(nsub-1)
    print 'rho(%d)=%g' % (ind1, sp.mean(t))

sp.savez('fcon1000_null_nsub50.npz', rho_null=rho_null, vrest_subs=vrest_subs)


# %%
# Read Candidate Subject to be tested against normals


rho_sub = sp.zeros([nsub, nVert])
for ind1 in range(nsub):
    vrest = vrest_subs[:, :, ind1]
    vrest, _, _ = normalizeData(vrest)
    vrest, Rot = brainSync(X=vsub, Y=vrest)
    t = sp.sum(vrest*vsub, axis=0)
#    print('rho(%d)=%g' % (ind1, sp.mean(t)), end=' ')
    print 'rho(%d)=%g' % (ind1, sp.mean(t))

    rho_sub[ind1, :] = t
    
    
# %%
# Hypothesis test

rho_null1 = sp.mean(rho_null, axis=1)
rho_sub1 = sp.mean(rho_sub, axis=0)
pval = sp.mean(rho_sub1 > rho_null1, axis=0)
r, corrPval,_,_ = multipletests(pvals=pval, alpha=0.05, method='fdr_bh')


sl = readdfs(os.path.join(BFPPATH, 'bci32kleft.dfs'))
sl.attributes = corrPval[:sl.vertices.shape[0]]
sl = patch_color_attrib(sl, clim=[0, 1])

sr = readdfs(os.path.join(BFPPATH, 'bci32kright.dfs'))
sr.attributes = corrPval[sl.vertices.shape[0]:2*sl.vertices.shape[0]]
sr = patch_color_attrib(sr, clim=[0, 1])

writedfs('right_pval_sn7602.dfs',sr);
writedfs('left_pval_sn7602.dfs',sl);

view_patch_vtk(sl, azimuth=90, elevation=180, roll=90, show=1)
view_patch_vtk(sl, azimuth=-90, elevation=180, roll=-90, show=1)
