% new full load

% load the options file (as opt here)
%load('/export/home/mike/matlab/batch/err_FPO92/fmri_err_cnth1s_FPO92_opt.mat')
%load('/export/home/mike/matlab/batch/err_FPO92_ds/fmri_err_cnth1s_FPO92_ds_opt.mat')
%load('/export/home/mike/matlab/batch/err_FPO82/fmri_err_cnth1s_FPO82_opt.mat') %this one used "pure coherence"
%load('/export/home/mike/matlab/batch/err_FPO82_freqds/fmri_err_cnth1s_FPO82_freqds_opt.mat')
load('/export/home/mike/matlab/batch/err_ctls_FPO82/err_ctl_cnth1s_ctls_FPO82_opt.mat')
%load('/export/home/mike/matlab/batch/err_FPO82_freqds_varcycle/fmri_err_cnth1s_FPO82_freqds_varcycle_opt.mat')
%load('/export/home/mike/matlab/batch/fmri_ern_beh.mat')
%load('/export/home/mike/matlab/batch/err_long/fmri_err_cnth1s_long_opt.mat')

%if we're using a legacy calculation (param_struct)
if exist('param_struct','var')
    outpath = '/active_projects/mike/fmri_phase4_ern_beh';
    coordsfile = '/export/home/mike/matlab/origin/coords/61chans_ns.mat';
    opt=fix_opt(param_struct,outpath,coordsfile);
    clear param_struct outpath coordsfile
end

%%

% set the demographics file if applicable
%demogsfile='/export/home/mike/matlab/origin/fmri/fmri_demogs_simple.mat';
demogsfile='/export/home/mike/matlab/csv/ph4_controls_over21_latestsesh.mat';
%demogsfile=[];

% set the factor for time-downsampling
time_ds=2;

% set the path for behavioral data if applicable
%behmat_path='/active_projects/mike/fmri_phase4_ern_beh';
behmat_path='/active_projects/mike/phase4_ern_beh';
%behmat_path=[];

%%

% import
%[mat_list, imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata, behdata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path);
%[mat_list, imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata, behdata, wave_evknormdata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path);
% valid types are
% n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata, behdata, wave_evknormdata
datatypes={'n_trials','erp','wave_evknorm','wave_tot','coh'};
[mat_list, imp, n_trials_all, erpdata, ~, wave_totdata, cohdata, ~, behdata, wave_evknormdata] = ...
    coh_import(opt, demogsfile, time_ds, behmat_path, datatypes);
%[mat_list, imp, n_trials_all, erpdata] = ...
%    coh_import(opt, demogsfile, time_ds);

%%

% scales
[scl,chan_locs] = coh_scales(opt,imp);

% rejection & demographics
trials_necessary=40;
%custom_rej=[7 10 14 23 28 41 52 54 57];
custom_rej=[];
[scl, s_inds_g, s_demogs] = coh_demog(mat_list, imp, scl, n_trials_all, demogsfile, ...
    trials_necessary,custom_rej);

% plotting parameters
[scl, pp] = coh_plotparams(imp,scl);

%%

% filter / baseline / pick peaks in ERP
lp_cutoff = 16;
[erpdata,peakmat] = coh_pickERPpeaks(imp, scl, s_inds_g, erpdata, [true true], lp_cutoff);

%%

clear datatypes demogsfile time_ds behmat_path lp_cutoff custom_rej trials_necessary