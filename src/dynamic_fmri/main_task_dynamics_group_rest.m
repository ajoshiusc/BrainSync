%||AUM||
%||Shree Ganeshaya Namaha||

clc;clear all;close all;
addpath(genpath('cifti_toolbox'));
addpath(genpath('myfuncs'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));

l=dir('/deneb_disk/HCP');
task='LANGUAGE';
subname='100307';
sub1zscore = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/tfMRI_%s/tfMRI_%s_hp200_s4_level2.feat/%s_tfMRI_%s_level2_hp200_s4.dscalar.nii',subname,task,task,subname,task);
sub1=ft_read_cifti(sub1zscore);
f1=fieldnames(sub1);

%zscore=sub1.x100307_tfmri_motor_level2_lh_hp200_s4
l=dir('/deneb_disk/HCP');
sub1zscore = ['/deneb_disk/HCP/135932/MNINonLinear/Results/tfMRI_',task,'/tfMRI_',task,'_hp200_s4_level2.feat/135932_tfMRI_',task,'_level2_hp200_s4.dscalar.nii'];
sub1=ft_read_cifti(sub1zscore);
imprv=[];
f=fieldnames(sub1);
avg_t1=cell(length(f),1);avg_t2=cell(length(f),1);avg_t2w=cell(length(f),1);
numVert=32492;
%load('/big_disk/ajoshi/with_andrew/reference/100307.reduce3.operators.mat')


%% Compute BrainSync Averaged task fmri
ref = '196750'; 
rfname = sprintf('/data_disk/HCP_All/%s/MNINonLinear/Results/tfMRI_%s_LR/tfMRI_%s_LR_Atlas.dtseries.nii',ref,task,task);
refdata=ft_read_cifti(rfname);
%rfname = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii',ref);
%refdata=ft_read_cifti(rfname);
refdata=normalizeData(refdata.dtseries');
tsum=0;nsub=0;
for subid=1:length(l)
        subname=l(subid).name;        
        subtdata = sprintf('/data_disk/HCP_All/%s/MNINonLinear/Results/tfMRI_%s_LR/tfMRI_%s_LR_Atlas.dtseries.nii',subname,task,task);
        
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
save(sprintf('tavg_%s_rest.mat',task),'dtseries');
%% Compute BrainSync Averaged rfmri
ref = '196750'; 
rfname = sprintf('/data_disk/HCP_All/%s/MNINonLinear/Results/tfMRI_%s_LR/tfMRI_%s_LR_Atlas.dtseries.nii',ref,task,task);
refdata=ft_read_cifti(rfname);
% rfname = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii',ref);
% refdata=ft_read_cifti(rfname);
refdata=normalizeData(refdata.dtseries');
rsum=0;nsub=0;
for subid=1:length(l)
        subname=l(subid).name;        
        subtdata = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii',subname);
        
        if ~isfile(subtdata)
            continue;
        end
        sub=ft_read_cifti(subtdata);
        rdata = sub.dtseries(:,1:316);
        rdata = normalizeData(rdata');
        rdata = brainSync(refdata,rdata);
        rsum = rsum+rdata;
        
        nsub=nsub+1
end
ravg=rsum/nsub;
save('ravg_nonorm.mat','ravg');

ravg=normalizeData(ravg);
save('ravg.mat','ravg');

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
    dd = brainSync(xin, D);
    dd = dd.*nrm;
    tskFitted(Cind+i-1, :) = dd(Cind, :);
    fprintf('%d,',i);
end
%[tskFitted,~, ~] = normalizeData(tskFitted);
save('tskFitted_rest.mat','tskFitted');

%%
zsum=0;
nsub=0;
for cind=7;%7%19 %1:length(f)
    if ~strcmp(f{cind}(end-2:end),'_s4')
        continue;
    end
    subind=1;
    for subid=1:length(l)
        subname=l(subid).name;
        
        sub2zscore = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/tfMRI_%s/tfMRI_%s_hp200_s4_level2.feat/%s_tfMRI_%s_level2_hp200_s4.dscalar.nii',subname,task,task,subname,task);
        
        if ~isfile(sub2zscore)
            continue;
        end
        sub2=ft_read_cifti(sub2zscore);
        zsum=zsum+sub2.([f{cind}(1),subname,f{cind}(8:end)]);
        nsub=nsub+1
    end
end
    
zavg=(zsum/nsub);
save('avg_zscore_rest.mat','zavg');
%% 
nVertHiRes=32492;
p_dir_ref = '/big_disk/ajoshi/HCP_data/';
ref = '196750'; 
%load('/big_disk/ajoshi/with_andrew/reference/100307.reduce3.operators.mat');

lsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.very_smooth.left.dfs']));
rsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.very_smooth.right.dfs']));

zl=zavg(1:nVertHiRes);%[ind_left]);
zr=zavg(nVertHiRes+1:2*nVertHiRes);%[ind_right]);
zscr=[zl;zr];


%% Show the average z scores
figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',1.0*(zl>5),'edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;

figure;
patch('faces',rsurf.faces,'vertices',rsurf.vertices,'facevertexcdata',1.0*(zr>5),'edgecolor','none','facecolor','interp');axis equal;axis off;view(90,0);
camlight; axis equal; axis off;material dull; 

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
save('tavg_rest.mat','tavg');
save('diffafter_rest.mat','dtseries');

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
save('diffafter_smooth_rest.mat','dtseries');

task_t=sqrt(sum(dtseries(:,tarea).^2,2));
% figure;
% plot(smooth(smooth(smooth(task_t))));

figure;
plot((task_t));


%%
