% new full load

% load the options file (as opt here)
load('/export/home/mike/matlab/batch/err_FPO92/fmri_err_cnth1s_FPO92_opt.mat')

% set the demographics file
demogsfile='/export/home/mike/matlab/origin/fmri/fmri_demogs_simple.mat';

% import
[mat_list, imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata] = ...
    coh_import_stage(opt, demogsfile);

% scales
[scl,chan_locs] = coh_scales_stage(opt,imp);

% rejection & demographics
[scl, s_inds_g, s_demogs] = coh_demog_stage(mat_list, imp, scl, n_trials_all, demogsfile);

% plotting parameters
pp = coh_plotparams_stage(imp,scl);

% filter / baseline / pick peaks in ERP
[erpdata,peakmat] = coh_pickERPpeaks_stage (imp, scl, s_inds_g, erpdata);