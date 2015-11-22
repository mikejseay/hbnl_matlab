% new full load

% load the options file (as opt here)
%load('/export/home/mike/matlab/batch/err_FPO92/fmri_err_cnth1s_FPO92_opt.mat')
%load('/export/home/mike/matlab/batch/err_FPO92_ds/fmri_err_cnth1s_FPO92_ds_opt.mat')
load('/export/home/mike/matlab/batch/err_FPO82/fmri_err_cnth1s_FPO82_opt.mat')
%load('/export/home/mike/matlab/batch/fmri_ern_beh.mat')

%if we're using a legacy calculation (param_struct)
if exist('param_struct','var')
    outpath = '/active_projects/mike/fmri_phase4_ern_beh';
    coordsfile = '/export/home/mike/matlab/origin/coords/61chans_ns.mat';
    opt=fix_opt(param_struct,outpath,coordsfile);
    clear param_struct outpath coordsfile
end

% set the demographics file if applicable
demogsfile='/export/home/mike/matlab/origin/fmri/fmri_demogs_simple.mat';

% set the factor for time-downsampling
time_ds=2;

% set the path for behavioral data if applicable
behmat_path='/active_projects/mike/fmri_phase4_ern_beh';

% import
%[mat_list, imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata, behdata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path);
[mat_list, imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata, behdata, wave_evknormdata] = ...
    coh_import(opt, demogsfile, time_ds, behmat_path);

% scales
[scl,chan_locs] = coh_scales(opt,imp);

% rejection & demographics
trials_necessary=15;
custom_rej=[7 10 14 23 28 41 52 57];
[scl, s_inds_g, s_demogs] = coh_demog(mat_list, imp, scl, n_trials_all, demogsfile, ...
    trials_necessary,custom_rej);

% plotting parameters
[scl, pp] = coh_plotparams(imp,scl);

% filter / baseline / pick peaks in ERP
lp_cutoff = 24;
[erpdata,peakmat] = coh_pickERPpeaks(imp, scl, s_inds_g, erpdata, lp_cutoff);

clear demogsfile time_ds behmat_path lp_cutoff custom_rej trials_necessary