# AUM
# Shree Ganeshaya Namaha
from dfsio import readdfs
from os.path import join
import nilearn.image
import numpy as np
from brainsync import normalizeData, brainSync
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


BFPPATH = '/big_disk/ajoshi/coding_ground/bfp'
BrainSuitePath = '/home/ajoshi/BrainSuite17a/svreg'

NCMP = 31

surfObj = readdfs(join(BFPPATH, 'supp_data', 'bci32kleft.dfs'))
numVert = len(surfObj.vertices)

sub1 = '/big_disk/ajoshi/HCP100/HCP100/135932/MNINonLinear/Results/rfMRI_REST1\
_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii'
sub1 = nilearn.image.load_img(sub1)
X = sub1.get_data().T
Xorig = np.array(X)
X, _, _ = normalizeData(X.T)

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

'''
D1 = D[:, :100]
X1 = X[:, :100]
p1=plt.figure()
plt.imshow(np.dot((D1.T), D1), aspect='auto', clim=(0, 1))
plt.title('estim')


p2=plt.figure()
plt.imshow(np.dot((X1.T), X1), aspect='auto', clim=(0, 1))
plt.title('orig')
# D is the representative data that produces the network
#DR = p.fit_transform(XR.T)
'''
nT = X.shape[0]

Xnew = np.zeros(X.shape)
Cind = np.ind((NCMP-1)/2+1)
for i in range(X.shape[0]-NCMP):
    xin = X[i:i+NCMP, :]
    xin, _, nrm = normalizeData(xin)
    dd, _ = brainSync(xin, D)
    dd = dd*nrm
    Xnew[Cind, :] = dd[Cind, :]
    print('%d/%d\n' % i, X.shape[0])

'''
% 
% %Xnew2=Xnew;
% %Xnew2(1:NCMP,:)=Xnew(1:NCMP,:)./(1:NCMP)';
% 
% %Xnew2(NCMP+1:end,:)=Xnew2(NCMP+1:end,:)/NCMP;
% c=sum(Xnew2.*X,1);
% cc=zeros(size(surfObj.vertices,1),1);
% cc(ind)=c;
% s=smooth_cortex_fast(surfObj,.1,600);
% 
% figure;
% patch('faces',surfObj.faces,'vertices',s.vertices,'facevertexcdata',cc,'edgecolor','none','facecolor','interp');
% axis equal;view(-90,0);camlight;material dull; axis off;caxis([0,1]);colormap jet;
% 
 '''
