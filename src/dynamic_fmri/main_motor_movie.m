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
zscore=sub1.x100307_tfmri_motor_level2_lh_hp200_s4
load('/big_disk/ajoshi/with_andrew/reference/100307.reduce3.operators.mat')

zl=zscore([ind_left]);
zr=zscore([ind_right]);
zscr=[zl;zr];
f1{7}


p_dir_ref = '/big_disk/ajoshi/HCP_data/'
ref = '196750' 
lsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.reduce3.very_smooth.left.dfs']));
rsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.reduce3.very_smooth.right.dfs']));


load('motor_diff_data_filt.mat')

figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',zl,'edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;

figure;
patch('faces',rsurf.faces,'vertices',rsurf.vertices,'facevertexcdata',zr,'edgecolor','none','facecolor','interp');axis equal;axis off;view(90,0);
camlight; axis equal; axis off;material dull;


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

