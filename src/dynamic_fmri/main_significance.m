%||AUM||
%||Shree Ganeshaya Namaha||

clc;clear all;close all;
addpath(genpath('cifti_toolbox'));
addpath(genpath('myfuncs'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
fmridatfile='/home/ajoshi/coding_ground/brainsync/src/dynamic_fmri/tavg.mat';
load(fmridatfile);
fmridatfile='/home/ajoshi/coding_ground/brainsync/src/dynamic_fmri/ravg.mat';
load(fmridatfile);
NCMP=21;
[~,D]=pca(ravg');
D=D(:,1:NCMP);
D=normalizeData(D');


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

rFitted = zeros(size(ravg));
Cind = (NCMP-1)/2+1;
for i = 1:size(ravg,1)-NCMP
    xin = ravg(i:i+NCMP-1, :);
    [xin, ~, nrm] = normalizeData(xin);
    dd = brainSync(xin, D);
    dd = dd.*nrm;
    rFitted(Cind+i-1, :) = dd(Cind, :);
    fprintf('%d,',i);
end

pval=0*ravg;
r1=ravg-rFitted;
t1=tavg-tskFitted;

nT=size(tavg,1);
for t=1:nT
%     parfor jj=1:size(tavg,2)
%         pval(t,jj)=sum(r1(:,jj)>t1(t,jj))/nT;
%     end
    parfor jj=1:size(tavg,2)
        pval(t,jj)=tcdf((mean(r1(:,jj))-t1(t,jj))/std(r1(:,jj)),1); % This is from 1 sample t-test formula from wikipedia https://en.wikipedia.org/wiki/Student's_t-test
    end
    t
end

p_dir_ref = '/big_disk/ajoshi/HCP_data/';
ref = '196750';


lsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.very_smooth.left.dfs']));
rsurf = readdfs(fullfile(p_dir_ref,'reference',[ref,'.aparc.a2009s.32k_fs.very_smooth.right.dfs']));
t=101;
nV=length(lsurf.vertices);
blk=79:105;% present story
%blk=[79:105,
blk=157:183;% present story
blk=268:299;%present story
%blk=137:145;% present math 
%blk=194:202;%present math
blk=175:177;
figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',1.0-(mean(pval(blk,1:nV))'),'edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;colormap jet

figure;
patch('faces',rsurf.faces,'vertices',rsurf.vertices,'facevertexcdata',1.0-(mean(pval(blk,1+nV:2*nV))'),'edgecolor','none','facecolor','interp');axis equal;axis off;view(90,0);
camlight; axis equal; axis off;material dull; colormap jet


save pval pval




%[tskFitted,~, ~] = normalizeData(tskFitted);
