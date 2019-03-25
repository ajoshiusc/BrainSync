%||AUM||
%||Shree Ganeshaya Namaha||
clc;clear all;close all;

NCMP=5;
% generate templet D from resting state
W=([1:(NCMP+1)/2,(NCMP-1)/2:-1:1].^2)';

load('/ImagePTE1/ajoshi/For_Anand_MICCAI/IdxNaN.mat')

subid = '100307';

fmri_rest_file=fullfile('/ImagePTE1/ajoshi/For_Anand_MICCAI/',subid,'rfMRI_REST1_LR.mat');

fmri_task_file=fullfile('/ImagePTE1/ajoshi/For_Anand_MICCAI/',subid,'tfMRI_LANGUAGE_LR.mat');

a=load(fmri_rest_file);
fmri_rest_origL = a.dataL;
fmri_restL=fmri_rest_origL(~idxNaNL,:);
fmri_restL=normalizeData(fmri_restL')';

fmri_rest_origR = a.dataR;
fmri_restR=fmri_rest_origR(~idxNaNR,:);
fmri_restR=normalizeData(fmri_restR')';

t=load(fmri_task_file);
fmri_task_origL = t.dataL;
fmri_taskL=fmri_task_origL(~idxNaNL,:);
fmri_taskL=normalizeData(fmri_taskL');

fmri_task_origR = t.dataR;
fmri_taskR=fmri_task_origR(~idxNaNR,:);
fmri_taskR=normalizeData(fmri_taskR');

dataL=0*fmri_task_origL;
dataR=0*fmri_task_origR;


[~,DL]=pca(fmri_restL);
DL=DL(:,1:NCMP);
DL=normalizeData(DL');

[~,DR]=pca(fmri_restR);
DR=DR(:,1:NCMP);
DR=normalizeData(DR');


fmri_task_fittedL = zeros(size(fmri_taskL));
fmri_task_fittedR = zeros(size(fmri_taskR));

Cind = (NCMP-1)/2+1;
for i = 1:size(fmri_taskL,1)-NCMP
    xinL = fmri_taskL(i:i+NCMP-1, :);
%    [xinL, ~, nrmL] = normalizeData(xinL);
    ddL = DbrainSync_wt(xinL, DL, W);
%    ddL = ddL.*nrmL;
    fmri_task_fittedL(Cind+i-1, :) = ddL(Cind, :);
    
    xinR = fmri_taskR(i:i+NCMP-1, :);
%    [xinR, ~, nrmR] = normalizeData(xinR);
    ddR = DbrainSync_wt(xinR, DR, W);
%    ddR = ddR.*nrmR;
    fmri_task_fittedR(Cind+i-1, :) = ddR(Cind, :);
    
    fprintf('%d,',i);
end
dataL(~idxNaNL,:)=fmri_task_fittedL';
dataR(~idxNaNR,:)=fmri_task_fittedR';

save('fmri_task_fitted_rest_wt.mat','dataL', 'dataR', '-v7.3');

