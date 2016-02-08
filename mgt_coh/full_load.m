% new full load

% load the options file (as opt here)
%load('/export/home/mike/matlab/batch/err_FPO92/fmri_err_cnth1s_FPO92_opt.mat')
%load('/export/home/mike/matlab/batch/err_FPO92_ds/fmri_err_cnth1s_FPO92_ds_opt.mat')
%load('/export/home/mike/matlab/batch/err_FPO82/fmri_err_cnth1s_FPO82_opt.mat') %this one used "pure coherence"
%load('/export/home/mike/matlab/batch/err_FPO82_freqds/fmri_err_cnth1s_FPO82_freqds_opt.mat')
%load('/export/home/mike/matlab/batch/err_ctls_FPO82/err_ctl_cnth1s_ctls_FPO82_opt.mat')
%load('/export/home/mike/matlab/batch/err_csub_FPO90/err_csubs_csub_FPO90_opt.mat')
load('/export/home/mike/matlab/batch/fmri_err_FPO90/fmri_err_cnth1s_fmri_err_FPO90_opt.mat') %fmri
%load /active_projects/Kam/MATLAB_KAM/MikeCoh/opt.mat %kam
%load('/export/home/mike/matlab/batch/err_csub_325/err_csubs_csub_325_opt.mat')
%load('/export/home/mike/matlab/batch/ern_csub_long/ern_csubs_FPO90_long_opt.mat')
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
% should have UniqueFileString, POP, and EEG_Age

%demogsfile='/export/home/mike/matlab/origin/fmri/fmri_demogs_simple.mat';
%demogsfile='/export/home/mike/matlab/csv/ph4_controls_over21_latestsesh.mat';
%demogsfile='/export/home/mike/matlab/batch/err_csub_FPO90/err_csubs_csub_FPO90_demog.mat';
%demogsfile='/export/home/mike/matlab/csv/kam_err_ph4_over21_ingroup.mat';
demogsfile='/export/home/mike/matlab/batch/fmri_err_FPO90/fmri_err_FPO90_demogs.mat';
%demogsfile=[];

% set the factor for time-downsampling
time_ds=4;

% set the path for behavioral data if applicable
%behmat_path='/active_projects/mike/fmri_phase4_ern_beh';
%behmat_path='/active_projects/mike/phase4_ern_beh';
%behmat_path='/active_projects/mike/ern_csub_beh';
behmat_path=[];

%%

% import
%[mat_list, imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata, behdata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path);
%[mat_list, imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata, behdata, wave_evknormdata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path);

% valid types are
% n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata, behdata, wave_evknormdata
%datatypes={'n_trials','erp','wave_evknorm','wave_tot','coh'};
%datatypes={'n_trials','erp','wave_totpow'};
%datatypes={'n_trials','erp','wave_evknorm','coh'};
datatypes={'n_trials','erp','wave_evknorm','coh'};

%[mat_list, imp, n_trials_all, behdata, erpdata, ~, wave_evknormdata, cohdata, phidata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path, datatypes);
% LONG ERP
%[mat_list, imp, n_trials_all, behdata, erpdata, wave_totpowdata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path, datatypes);
% REGULAR
%[mat_list, imp, n_trials_all, behdata, erpdata, ~, wave_evknormdata, cohdata, phidata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path, datatypes);
% KAM
[mat_list, imp, n_trials_all, ~, erpdata, ~, wave_evknormdata, cohdata, ~, ~] = ...
    coh_import(opt, demogsfile, time_ds, behmat_path, datatypes);


%%

% scales
[scl,chan_locs] = coh_scales(opt,imp);

% rejection & demographics
trials_necessary=40; %15 for ERP/power, 40 for ITC/ISPC
%custom_rej=[7 10 14 23 28 41 52 54 57];
%custom_rej=[12,28,38,39,41,43,61,63,65,68,82,86,87,89,97,99,101,108,115,117,119,124,126,129,133,134,135,146,157,159,160,161];
%custom_rej=[1,5,6,8,14,22,23,26,27,33,40,41,45,47,54,55,58,60,64,69,78,82,93,102,103,110,115,118,127,147,158,167,168,169,175,182,189,194];
%custom_rej= [230,242,540,545,1044,1345,1457,1464,1468,1483,1492,1825,1826,1839,2030,2034];
custom_rej=[];
%age_range=[21 45];
age_range=[];
[scl, s_inds_g, s_demogs] = coh_demog(mat_list, imp, scl, n_trials_all, demogsfile, ...
    trials_necessary, age_range, custom_rej);

% plotting parameters
cond_diff = {2, 1};
[scl, pp] = coh_plotparams(imp, scl, cond_diff);

%%

% filter / baseline / pick peaks in ERP
lp_cutoff = 16;
[erpdata,peakmat] = coh_pickERPpeaks(imp, scl, s_inds_g, erpdata, [true true], lp_cutoff);

%%

clear datatypes demogsfile time_ds behmat_path lp_cutoff custom_rej ...
    trials_necessary age_range cond_diff