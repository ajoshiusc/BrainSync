%||AUM||
%||Shree Ganeshaya Namaha||

clc;clear all;close all;
addpath(genpath('cifti_toolbox'));
addpath(genpath('myfuncs'));
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'));
fmridatfile='/home/ajoshi/coding_ground/brainsync/src/dynamic_fmri/tavg_LANGUAGE_rest_wt.mat';

load(fmridatfile);tavg=dtseries;
fmridatfile='/home/ajoshi/coding_ground/brainsync/src/dynamic_fmri/ravg_wt.mat';
load(fmridatfile);
NCMP=21;
W=([1:(NCMP+1)/2,(NCMP-1)/2:-1:1].^2)';
%W=ones(size(W));
[~,D]=pca(ravg');
D=D(:,1:NCMP);
D=normalizeData(D');


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

rFitted = zeros(size(ravg));
Cind = (NCMP-1)/2+1;
for i = 1:size(ravg,1)-NCMP
    xin = ravg(i:i+NCMP-1, :);
    [xin, ~, nrm] = normalizeData(xin);
    dd = brainSync_wt(xin, D, W);
    dd = dd.*nrm;
    rFitted(Cind+i-1, :) = dd(Cind, :);
    fprintf('%d,',i);
end

pval=0*ravg;
pval_rest=0*ravg;

r1=ravg-rFitted;
t1=tavg-tskFitted;

nT=size(tavg,1);
for t=1:nT
%     parfor jj=1:size(tavg,2)
%         pval(t,jj)=sum(r1(:,jj)>t1(t,jj))/nT;
%     end
    parfor jj=1:size(tavg,2)
        pval(t,jj)=tcdf(mean(r1(:,jj)-t1(t,jj))/std(r1(:,jj)),1); % This is from 1 sample t-test formula from wikipedia https://en.wikipedia.org/wiki/Student's_t-test
        pval_rest(t,jj)=tcdf(mean(r1(:,jj)-r1(t,jj))/std(r1(:,jj)),1); % This is from 1 sample t-test formula from wikipedia https://en.wikipedia.org/wiki/Student's_t-test
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
%blk=283:287;
%blk=185:187; % response story
figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',1.0-(mean(pval(blk,1:nV))'),'edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;colormap jet;axis tight;

figure;
patch('faces',rsurf.faces,'vertices',rsurf.vertices,'facevertexcdata',1.0-(mean(pval(blk,1+nV:2*nV))'),'edgecolor','none','facecolor','interp');axis equal;axis off;view(90,0);
camlight; axis equal; axis off;material dull; colormap jet;axis tight;

figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',1.0-(mean(pval_rest(blk,1:nV))'),'edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;colormap jet;axis tight;

figure;
patch('faces',rsurf.faces,'vertices',rsurf.vertices,'facevertexcdata',1.0-(mean(pval_rest(blk,1+nV:2*nV))'),'edgecolor','none','facecolor','interp');axis equal;axis off;view(90,0);
camlight; axis equal; axis off;material dull; colormap jet;axis tight;
%%%%
figure;
patch('faces',lsurf.faces,'vertices',lsurf.vertices,'facevertexcdata',mean(tskFitted(blk,1:nV),1)','edgecolor','none','facecolor','interp');axis equal;axis off;view(-90,0);
camlight; axis equal; axis off;material dull;colormap jet;axis tight;



save pval_wt pval pval_rest

pval_sm=0*pval;
parfor jj=1:size(pval,2)
    pval_sm(:,jj)=smooth(pval(:,jj));
end

save pval_sm_wt pval_sm


%[tskFitted,~, ~] = normalizeData(tskFitted);
