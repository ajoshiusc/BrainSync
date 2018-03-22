%||AUM||
%||Shree Ganeshaya Namaha||

clc;clear all;close all;
addpath(genpath('cifti_toolbox'));
addpath(genpath('myfuncs'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));

l=dir('/deneb_disk/HCP');
numVert=32492;
%load('/big_disk/ajoshi/with_andrew/reference/100307.reduce3.operators.mat')

%% Compute BrainSync Averaged rfmri
ref = '196750';
rfname = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii',ref);
refdata=ft_read_cifti(rfname);
refdata=normalizeData(refdata.dtseries');
rsum=0;nsub=0;
for subid=1:length(l)
    subname=l(subid).name;
    subtdata = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii',subname);
    
    if ~isfile(subtdata)
        continue;
    end
    sub=ft_read_cifti(subtdata);
    rdata = sub.dtseries;
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
[~,D]=pca(ravg');
D=D(:,1:NCMP);
D=normalizeData(D');

%% Compute Dynamics of task data

rstFitted = zeros(size(ravg));
Cind = (NCMP-1)/2+1;
diff=zeros(size(refdata,1),size(refdata,2));
nsub=0;
for subid=1:length(l)
    subname=l(subid).name;
    subtdata = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii',subname);
    
    if ~isfile(subtdata)
        continue;
    end
    sub=ft_read_cifti(subtdata);
    rdata = sub.dtseries;
    rdata = normalizeData(rdata');
    rdata = brainSync(refdata,rdata);
    
    
    for i = 1:size(ravg,1)-NCMP
        xin = rdata(i:i+NCMP-1, :);
        [xin, ~, nrm] = normalizeData(xin);
        dd = brainSync(xin, D);
        dd = dd.*nrm;
        rstFitted(Cind+i-1, :) = dd(Cind, :);
        diff(Cind+i-1,:)=diff(Cind+i-1,:)+(rdata(Cind+i-1, :)-rstFitted(Cind+i-1, :)).^2;
        
    end
    nsub=nsub+1;
    fprintf('%d,',nsub);
    
end
%save('rstFitted.mat','rstFitted');
save('diffrest.mat', 'diff');
save('diffrest2.mat');

%%
nVertHiRes=32492;
diffn=sqrt(sum(diff.^2,1));
difft=sqrt(sum(diff.^2,2));

p_dir_ref = '/big_disk/ajoshi/HCP_data/';
ref = '196750';
%load('/big_disk/ajoshi/with_andrew/reference/100307.reduce3.operators.mat');

lsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.very_smooth.left.dfs']));
rsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.very_smooth.right.dfs']));



h=figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',diffn(1:length(lsurf.vertices))','edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
axis equal; axis off;material dull;colormap jet;caxis([.5,2.5]);camlight;
saveas(h,'fig1.fig');
h2=figure;
patch('faces',rsurf.faces,'vertices',rsurf.vertices,'facevertexcdata',diffn(1+length(lsurf.vertices):2*length(lsurf.vertices))','edgecolor','none','facecolor','interp');axis equal;axis off;view(90,0);
axis equal; axis off;material dull;colormap jet;caxis([1.5,2.5]);camlight;
saveas(h2,'fig2.fig');
