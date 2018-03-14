%||AUM||
%||Shree Ganeshaya Namaha||

clc;clear all;close all;
addpath(genpath('cifti_toolbox'));
addpath(genpath('myfuncs'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));

l=dir('/deneb_disk/HCP');
task='MOTOR';
subname='100307';
sub1zscore = sprintf('/deneb_disk/HCP/%s/MNINonLinear/Results/tfMRI_%s/tfMRI_%s_hp200_s4_level2.feat/%s_tfMRI_%s_level2_hp200_s4.dscalar.nii',subname,task,task,subname,task);
sub1=ft_read_cifti(sub1zscore);
f1=fieldnames(sub1);

%zscore=sub1.x100307_tfmri_motor_level2_lh_hp200_s4
l=dir('/deneb_disk/HCP');
task='MOTOR';'WM';'LANGUAGE';
sub1zscore = ['/deneb_disk/HCP/135932/MNINonLinear/Results/tfMRI_',task,'/tfMRI_',task,'_hp200_s4_level2.feat/135932_tfMRI_',task,'_level2_hp200_s4.dscalar.nii'];
sub1=ft_read_cifti(sub1zscore);
imprv=[];
f=fieldnames(sub1);
avg_t1=cell(length(f),1);avg_t2=cell(length(f),1);avg_t2w=cell(length(f),1);
numVert=32492;
load('/big_disk/ajoshi/with_andrew/reference/100307.reduce3.operators.mat')
zsum=0;
nsub=0;
for cind=4;%7%19 %1:length(f)
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

%% 
nVertHiRes=32492;
p_dir_ref = '/big_disk/ajoshi/HCP_data/'
ref = '196750' 
load('/big_disk/ajoshi/with_andrew/reference/100307.reduce3.operators.mat')

lsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.very_smooth.left.dfs']));
rsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.very_smooth.right.dfs']));

zl=zavg(1:nVertHiRes);%[ind_left]);
zr=zavg(nVertHiRes+1:2*nVertHiRes);%[ind_right]);
zscr=[zl;zr];


%% Show the average z scores
figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',zl,'edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;

figure;
patch('faces',rsurf.faces,'vertices',rsurf.vertices,'facevertexcdata',zr,'edgecolor','none','facecolor','interp');axis equal;axis off;view(90,0);
camlight; axis equal; axis off;material dull; 






%% 
diffrt=sum(diffafter.^2,1);
figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',diffrt(1:length(lsurf.vertices))','edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;

figure;
patch('faces',rsurf.faces,'vertices',rsurf.vertices,'facevertexcdata',diffrt(1+length(rsurf.vertices):end)','edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;


tarea=find((zscr>5));

diffafter(:,tarea)
tongue_t=sqrt(sum(diffafter(:,tarea).^2,2));
figure;
plot(smooth(smooth(tongue_t)));

