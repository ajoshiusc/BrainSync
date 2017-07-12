clc;clear all;close all;
addpath(genpath('/big_disk/ajoshi/coding_ground/svreg-matlab/src'))
addpath(genpath('/big_disk/ajoshi/coding_ground/svreg-matlab/3rdParty'))
addpath(genpath('/big_disk/ajoshi/coding_ground/brainreg'))
h_var=5;


l=dir('/big_disk/ajoshi/fcon_1000/Beijing/sub*');
p_dir = '/big_disk/ajoshi/with_andrew';
ref_dir = fullfile(p_dir, 'reference');
r_factor = 3;
N = 11;
ref = '100307';
conn = load(fullfile(ref_dir, [ref '.reduce' num2str(r_factor) '.vertex_conn_' num2str(N) 'N.mat']));
for subno = 1:length(l)
    fname = l(subno).name;
    subno
    if ~exist(['/big_disk/ajoshi/fcon_1000/Beijing/',fname,'/fmrit_reduce3_v2.mat'],'file')        
        continue;
    end
    if exist(['/big_disk/ajoshi/fcon_1000/Beijing/',fname,'/fmri_tnlm_5_reduce3_v2.mat'],'file')
        continue;
    end
    
    load(['/big_disk/ajoshi/fcon_1000/Beijing/',fname,'/fmrit_reduce3_v2.mat']);
    
    
    func_left = tNLM_SpatioTemporalData(fmri_left, conn.v_conn_left, h_var);
    
    func_right = tNLM_SpatioTemporalData(fmri_right, conn.v_conn_right, h_var);
    
    save(['/big_disk/ajoshi/fcon_1000/Beijing/',fname,'/fmri_tnlm_5_reduce3_v2.mat'], 'func_left', 'func_right');
end

