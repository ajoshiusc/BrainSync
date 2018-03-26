%||AUM||
%||Shree Ganeshaya Namaha||

clc;clear all;close all;
addpath(genpath('cifti_toolbox'));
addpath(genpath('myfuncs'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
fmridatfile='/home/ajoshi/coding_ground/brainsync/src/dynamic_fmri/diffafter_rest.mat';
%fmridatfile='/home/ajoshi/coding_ground/brainsync/src/dynamic_fmri/tavg_rest.mat';

%fmridatfile='/home/ajoshi/coding_ground/brainsync/src/dynamic_fmri/tskFitted.mat';
clim=[-.1,.1];
cmap=bipolarcmapW(100,clim,'linear','br');
cname='diffafter_rest';


fl1=sprintf('%s_l1.png',cname);fl2=sprintf('%s_l2.png',cname);
fr1=sprintf('%s_r1.png',cname);fr2=sprintf('%s_r2.png',cname);

load(fmridatfile);

if exist('tavg','var')
    dtseries=tavg;
%     for jj=1:size(dtseries,2)
%         dtseries(:,jj)=smooth(dtseries(:,jj));
%     end
end

if exist('tskFitted','var')
    dtseries=tskFitted;
end


nVertHiRes=32492;
diffn=sqrt(sum(dtseries.^2,1));
difft=sqrt(sum(dtseries.^2,2));

p_dir_ref = '/big_disk/ajoshi/HCP_data/';
ref = '196750';
%load('/big_disk/ajoshi/with_andrew/reference/100307.reduce3.operators.mat');

lsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.very_smooth.left.dfs']));
rsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.very_smooth.right.dfs']));

lav=0*(mean(dtseries,1));
blk=79:105;% present story
blk=157:183;% present story
%blk=268:299;%present story
%blk=137:145;% present math 
%blk=194:202;%present math
%blk=150:154;
blk=176;
l=(sum((dtseries(blk,:)-lav),1));

l1=l(1:length(lsurf.vertices))';%smooth_surf_function(lsurf,l(1:length(lsurf.vertices))',.3,.3);
h=figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',l1,'edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;axis tight;colormap(cmap);axis tight;caxis(clim);
saveas(h,fl1);view(90,0);camlight;saveas(h,fl2);
%caxis([-2,2]);

l2=l((1+length(rsurf.vertices)):2*length(rsurf.vertices))';%smooth_surf_function(rsurf,l((1+length(rsurf.vertices)):2*length(rsurf.vertices))',.3,.3);
h=figure;
patch('faces',rsurf.faces,'vertices',rsurf.vertices,'facevertexcdata',l2,'edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
view(90,0);camlight; axis equal; axis off;material dull;axis tight;colormap(cmap);axis tight;caxis(clim);
saveas(h,fr1);view(-90,0);camlight;saveas(h,fr2);
%caxis([-2,2]);

