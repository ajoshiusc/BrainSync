# AUM
# Shree Ganeshaya Namaha
from dfsio import readdfs
from os.path import join
import nilearn.image
import numpy as np
from brainsync import normalizeData, brainSync
from sklearn.decomposition import PCA
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
#XL = X[:numVert, :]
#XR = X[1+numVert:2*numVert, :]
#indL = np.isfinite(np.sum(XL, axis=1))
#indR = np.isfinite(np.sum(XR, axis=1))

#XL = XL[indL, :]
#XR = XR[indR, :]
X, _, _ = normalizeData(X.T)
#XR, _, _ = normalizeData(XR.T)

p = PCA(n_components=NCMP)

D = p.fit_transform(X.T)

D = normalizeData(D.T)

#DR = p.fit_transform(XR.T)

'''
_,SL,VL=svd(XL);[~,SR,VR]=svd(XR);

% This is a representative data that produces the network.
DL=SL(1:NCMP,1:NCMP)*VL(:,1:NCMP)';DR=SR(1:NCMP,1:NCMP)*VR(:,1:NCMP)';

%normalize to unit norm?
%don't subtract the mean
clear SL SR VL VR sub1

nTL=size(XL,1);nTR=size(XR,1);

XnewL=zeros(size(XL,1),size(XL,2));XnewR=zeros(size(XR,1),size(XR,2));
Xnew2L=zeros(size(XL));Xnew2R=zeros(size(XR));
dL=sqrt(sum(DL.^2,1));dR=sqrt(sum(DR.^2,1));
D1L=DL./dL; D1R=DR./dR;

for i=1:size(XL,1)-NCMP
    xinL=XL(i:i+NCMP-1,:);xinR=XR(i:i+NCMP-1,:);
    dL=sqrt(sum(xinL.^2,1));dR=sqrt(sum(xinR.^2,1));
    xinL=xinL./dL;xinR=xinR./dR;
    ddL= brainSync(xinL,D1L);ddL=ddL.*dL; ddR= brainSync(xinR,D1R);ddR=ddR.*dR;
    Xnew2L(i+(NCMP-1)/2,:) =  ddL((NCMP-1)/2+1,:); Xnew2R(i+(NCMP-1)/2,:) =  ddR((NCMP-1)/2+1,:);
    fprintf('%d/%d\n',i,size(XL,1));
end
Xorig2=normalizeData(Xorig')';
iiL=find(indL);iiR=find(indR)+numVert;
XX=Xorig2;
dtseries=XX;
save orig_dyn.mat dtseries


XX(iiL,:)=Xnew2L';XX(iiR,:)=Xnew2R';
dtseries=XX;
save fitted_dyn.mat dtseries
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
