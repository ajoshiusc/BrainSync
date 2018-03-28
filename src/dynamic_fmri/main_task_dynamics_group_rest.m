%||AUM||
%||Shree Ganeshaya Namaha||

clc;clear all;close all;
addpath(genpath('cifti_toolbox'));
addpath(genpath('myfuncs'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
NCMP=21;

W=([1:(NCMP+1)/2,(NCMP-1)/2:-1:1].^2)';
l=dir('/deneb_disk/HCP');
task='LANGUAGE';
subname='100307';
sub1zscore = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/tfMRI_%s/tfMRI_%s_hp200_s4_level2.feat/%s_tfMRI_%s_level2_hp200_s4.dscalar.nii',subname,task,task,subname,task);
sub1=ft_read_cifti(sub1zscore);
f1=fieldnames(sub1);

l=dir('/deneb_disk/HCP');
numVert=32492;


%% Compute BrainSync Averaged task fmri
ref = '196750'; 
task = 'LANGUAGE';
rfname = sprintf('/data_disk/HCP_All/%s/MNINonLinear/Results/tfMRI_%s_LR/tfMRI_%s_LR_Atlas.dtseries.nii',ref,task,task);
r=ft_read_cifti(rfname);
%r=load(sprintf('/deneb_disk/HCP_filt_data/%s_tfMRI_%s_LR_Atlas.dtseries.filt.mat',ref,task));
refdata=r.dtseries;
refdata=normalizeData(refdata');

tsum=0;nsub=0;
for subid=1:length(l)
        subname=l(subid).name;        
        subtdata = sprintf('/data_disk/HCP_All/%s/MNINonLinear/Results/tfMRI_%s_LR/tfMRI_%s_LR_Atlas.dtseries.nii',subname,task,task);
%        subtdata = sprintf('/deneb_disk/HCP_filt_data/%s_tfMRI_%s_LR_Atlas.dtseries.filt.mat',subname,task);
        
        if ~isfile(subtdata)
            continue;
        end
        sub=ft_read_cifti(subtdata);
        tdata = sub.dtseries;
        tdata = normalizeData(tdata');
        tdata = brainSync(refdata,tdata);
        tsum = tsum+tdata;
        
        nsub=nsub+1
end
tavg=tsum/nsub;
tavg=normalizeData(tavg);dtseries=tavg;
save(sprintf('tavg_%s_rest_wt.mat',task),'dtseries');
%% Compute BrainSync Averaged rfmri
ref = '196750'; 
rfname = sprintf('/data_disk/HCP_All/%s/MNINonLinear/Results/tfMRI_%s_LR/tfMRI_%s_LR_Atlas.dtseries.nii',ref,task,task);
%rfname = sprintf('/deneb_disk/HCP_filt_data/%s_rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.filt.mat',ref);
%refdata=load(rfname);
refdata=ft_read_cifti(rfname);
refdata=refdata.dtseries(:,1:316);
refdata=normalizeData(refdata');
rsum=0;nsub=0;
for subid=1:length(l)
        subname=l(subid).name;        
        subtdata = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii',subname);
%        subtdata = sprintf('/deneb_disk/HCP_filt_data/%s_rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.filt.mat',subname);
        
        if ~isfile(subtdata)
            subtdata
            continue;
            
        end
        sub=ft_read_cifti(subtdata);
        %sub=load(subtdata);
        rdata = sub.dtseries(:,1:316);
        rdata = normalizeData(rdata');
        rdata = brainSync(refdata,rdata);
        rsum = rsum+rdata;
        
        nsub=nsub+1
end
ravg=rsum/nsub;
save('ravg_nonorm_wt.mat','ravg');

ravg=normalizeData(ravg);
save('ravg_wt.mat','ravg');

%% Compute the dynamics from resting data  
%% Alternatively, template can be built using task data, which is slightly worse
NCMP=21;
[~,D]=pca(tavg');
D=D(:,1:NCMP);
D=normalizeData(D');

%% Compute Dynamics of task data

tskFitted = zeros(size(tavg));
Cind = (NCMP-1)/2+1;
for i = 1:size(tavg,1)-NCMP
    xin = tavg(i:i+NCMP-1, :);
    [xin, ~, nrm] = normalizeData(xin);
    dd = brainSync_wt(xin, D, W);
    dd = dd.*nrm;
    tskFitted(Cind+i-1, :) = dd(Cind, :);
    fprintf('%d,',i);
end
%[tskFitted,~, ~] = normalizeData(tskFitted);
save('tskFitted_rest_wt.mat','tskFitted');

%% Compute and Plot difference

diffafter=tavg-tskFitted;

diffrt=sqrt(sum(diffafter.^2,1));

figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',diffrt(1:length(lsurf.vertices))','edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;

figure;
patch('faces',rsurf.faces,'vertices',rsurf.vertices,'facevertexcdata',diffrt((1+length(rsurf.vertices)):2*length(rsurf.vertices))','edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;

tarea=find((zavg>5));
dtseries=diffafter;
save('tavg_rest_wt.mat','tavg');
save('diffafter_rest_wt.mat','dtseries');

% diffafter(:,tarea);
% task_t=sqrt(sum(diffafter(:,tarea).^2,2));
% figure;
% plot(smooth(smooth(smooth(task_t))));
% 
% figure;
% plot(smooth(smooth((task_t))));
% 
dtseries=zeros(size(diffafter));
for jj=1:size(dtseries,2)
    dtseries(:,jj)=smooth(diffafter(:,jj));
end
save('diffafter_smooth_rest_wt.mat','dtseries');

task_t=sqrt(sum(dtseries(:,tarea).^2,2));
% figure;
% plot(smooth(smooth(smooth(task_t))));

figure;
plot((task_t));


%%
