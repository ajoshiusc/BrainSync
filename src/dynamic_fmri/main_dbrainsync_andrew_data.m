%||AUM||
%||Shree Ganeshaya Namaha||
clc;clear all;close all;

NCMP=21;
% generate templet D from resting state
W=([1:(NCMP+1)/2,(NCMP-1)/2:-1:1].^2)';

load('/ImagePTE1/ajoshi/For_Anand_MICCAI/IdxNaN.mat')

subid = '100307';

fmri_rest_file=fullfile('/ImagePTE1/ajoshi/For_Anand_MICCAI/',subid,'rfMRI_REST1_LR.mat');

fmri_task_file=fullfile('/ImagePTE1/ajoshi/For_Anand_MICCAI/',subid,'tfMRI_LANGUAGE_LR.mat');

a=load(fmri_rest_file);
fmri_rest_orig = a.dataL;
fmri_rest=fmri_rest_orig(~idxNaNL,:);
fmri_rest=normalizeData(fmri_rest')';


t=load(fmri_task_file);
fmri_task_orig = t.dataL;
fmri_task=fmri_task_orig(~idxNaNL,:);
fmri_task=normalizeData(fmri_task');

dataL=0*fmri_task_orig;


[~,D]=pca(fmri_rest);
D=D(:,1:NCMP);
D=normalizeData(D');


fmri_task_fitted = zeros(size(fmri_task));
Cind = (NCMP-1)/2+1;
for i = 1:size(fmri_task,1)-NCMP
    xin = fmri_task(i:i+NCMP-1, :);
    [xin, ~, nrm] = normalizeData(xin);
    dd = DbrainSync_wt(xin, D, W);
    dd = dd.*nrm;
    fmri_task_fitted(Cind+i-1, :) = dd(Cind, :);
    fprintf('%d,',i);
end
dataL(~idxNaNL,:)=fmri_task_fitted';
%[tskFitted,~, ~] = normalizeData(tskFitted);
save('fmri_task_fitted_rest_wt.mat','dataL');

