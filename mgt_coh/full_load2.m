% new full load

% load the options file (as opt here)
load('/export/home/mike/matlab/batch/err_FPO92/fmri_err_cnth1s_FPO92_opt.mat')

% set the demographics file
demogsfile='/export/home/mike/matlab/origin/fmri/fmri_demogs_simple.mat';

% import
[imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata] = ...
    coh_import_stage(opt, demogsfile);

% scales