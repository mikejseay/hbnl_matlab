
load fmri_ern_wavelet_cleanCSD_newconds2.mat
load('/export/home/mike/matlab/origin/mgt_coh/fmri_ern_h1_files.mat')

h1_struct = read_hdf1_dataV7(h1_list{2});
etable=h1_getbehav(h1_struct,param_struct.trial_init_type);
[param_struct.trial_mat,param_struct.case_label,respRT] = ...
behavinds_sog2(etable);
n_trials=sum(param_struct.trial_mat,1)