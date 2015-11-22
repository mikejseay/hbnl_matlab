%load h1 file list
load('/export/home/mike/matlab/origin/mgt_coh/fmri_err_h1_files.mat')
% params (basic)
load('/export/home/mike/matlab/batch/err_cleanCSD.mat')

data_file_type='hdf_binary';
info_file = [];
option = [false false false false];
vis = false;

output = erp_analysis3(h1_list{1}, data_file_type, info_file, ...
    param_struct, option, vis);