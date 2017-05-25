clc;clear all;close all;
addpath(genpath('/home/ajoshi/coding_ground/svreg/src'))
addpath(genpath('/home/ajoshi/coding_ground/svreg/3rdParty'))
addpath(genpath('/home/ajoshi/coding_ground/brainreg'))
h_var=0.15;


l=dir('/deneb_disk/studyforrest/sub*');
p_dir = '/big_disk/ajoshi/HCP_data';
ref_dir = fullfile(p_dir, 'reference');
r_factor = 3;
N = 11;
ref = '100307';
conn = load(fullfile(ref_dir, [ref '.reduce' num2str(r_factor) '.vertex_conn_' num2str(N) 'N.mat']));
sr=readdfs('/big_disk/ajoshi/HCP_data/reference/100307.aparc.a2009s.32k_fs.reduce3.very_smooth.right.dfs');
for subno = 1:length(l)
    fname = l(subno).name;
    subno
    if ~exist(['/deneb_disk/studyforrest/',fname,'/fmrit_reduce3_v2.mat'],'file')        
        continue;
    end
    if exist(['/deneb_disk/studyforrest/',fname,'/fmri_tnlm_0p15_reduce3_v2.mat'],'file')
        continue;
    end
    
    load(['/deneb_disk/studyforrest/',fname,'/fmrit_reduce3_v2.mat']);
    
    figure;patch('faces',sr.faces,'vertices',sr.vertices,'facevertexcdata',double(isnan(fmri_right(:,1))),'edgecolor','none','facecolor','interp');
    axis equal;axis off;camlight;material dull;
        
    func_left = tNLM_SpatioTemporalData(fmri_left, conn.v_conn_left, h_var);
    
    func_right = tNLM_SpatioTemporalData(fmri_right, conn.v_conn_right, h_var);
    
    save(['/deneb_disk/studyforrest/',fname,'/fmri_tnlm_0p15_reduce3_v2.mat'], 'func_left', 'func_right');
end




