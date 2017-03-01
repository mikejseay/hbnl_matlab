% new full load

% load the options file (as opt here)
%load('/export/home/mike/matlab/batch/err_FPO92/fmri_err_cnth1s_FPO92_opt.mat')
%load('/export/home/mike/matlab/batch/err_FPO92_ds/fmri_err_cnth1s_FPO92_ds_opt.mat')
%load('/export/home/mike/matlab/batch/err_FPO82/fmri_err_cnth1s_FPO82_opt.mat') %this one used "pure coherence"
%load('/export/home/mike/matlab/batch/err_FPO82_freqds/fmri_err_cnth1s_FPO82_freqds_opt.mat')
%load('/export/home/mike/matlab/batch/err_ctls_FPO82/err_ctl_cnth1s_ctls_FPO82_opt.mat')
%load('/export/home/mike/matlab/batch/err_csub_FPO90/err_csubs_csub_FPO90_opt.mat') %csubs % *****
%load('/export/home/mike/matlab/batch/fmri_err_FPO90/fmri_err_cnth1s_fmri_err_FPO90_opt.mat') %fmri
%load /active_projects/Kam/MATLAB_KAM/MikeCoh/opt.mat %kam
%load('/export/home/mike/matlab/batch/ashwini_err_fpo90/ashwini_err_h1sashwini_fpo90_opt.mat') %ashwini sample
%load('/export/home/mike/matlab/batch/asubs_err_fpo90/asubs_err_h1sasubs_fpo90_opt.mat')  % **** %asubs
%load('/export/home/mike/matlab/batch/err_csub_325/err_csubs_csub_325_opt.mat')
%load('/export/home/mike/matlab/batch/ern_csub_long/ern_csubs_FPO90_long_opt.mat')
%load('/export/home/mike/matlab/batch/err_FPO82_freqds_varcycle/fmri_err_cnth1s_FPO82_freqds_varcycle_opt.mat')
%load('/export/home/mike/matlab/batch/fmri_ern_beh.mat')
%load('/export/home/mike/matlab/batch/err_long/fmri_err_cnth1s_long_opt.mat')
%load('/export/home/mike/matlab/batch/fmri_err_fpo90_vcycle/fmri_err_cnth1sfpo90_vcycle_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_vp3_325/ac252_vp3_325_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_gng_fpo90/ac252_gng_fpo90_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_aod_fpo90/ac252_aod_fpo90_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_vp3_nonCSD/ac252_vp3nonCSD_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_cpt_nonCSD/ac252_cptnonCSD_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_vp3_nonCSD_mono/ac252_vp3_nonCSD_mono_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_cpt_nonCSD_mono/ac252_cpt_nonCSD_mono_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_ans_nonCSD_mono/ac252_ans_nonCSD_mono_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_stp_nonCSD_mono/ac252_stp_nonCSD_mono_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_vp3_noasr_mono/ac252_vp3ac252_vp3_noasr_mono_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_vp3_asr_mono/ac252_vp3ac252_vp3_asr_mono_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_vp3_asr_noburst/ac252_vp3ac252_vp3_asr_noburst_opt.mat')
%load('/export/home/mike/matlab/batch/nki_vp3_eog_lms/nki_vp3_cnth1s_eog_lms_opt.mat')
%load('/export/home/mike/matlab/batch/nki_vp3_noclean_noeog/nki_vp3_cnth1s_noclean_noeog_opt.mat')
%load('/export/home/mike/matlab/batch/nki_vp3_clean_noeog/nki_vp3_cnth1s_clean_noeog_opt.mat')
%load('/export/home/mike/matlab/batch/nki_vp3_noclean_noeog_simple/nki_vp3_cnth1s_noclean_noeog_simple_opt.mat')
%load('/export/home/mike/matlab/batch/nki_vp3_noclean_03_noeog_simple/nki_vp3_cnth1s_noclean_03_noeog_simple_opt.mat')
%load('/export/home/mike/matlab/batch/nki_err_nonCSD/nki_err_cnth1snonCSD_opt.mat')
%load('/export/home/mike/matlab/batch/nki_ern_nonCSD/nki_ern_cnth1snonCSD_opt.mat')
%load('/export/home/mike/matlab/batch/nki_err_CSD/nki_err_cnth1sCSD_opt.mat')
%load('/export/home/mike/matlab/batch/nki_vp3_bsseognew_both/nki_vp3_cnth1s_bsseognew_both_opt.mat')
%load('/export/home/mike/matlab/batch/nki_vp3_bsseog_both/nki_vp3_cnth1s_bsseog_both_opt.mat')
%load('/export/home/mike/matlab/batch/ac252_err_bssCSD/ac252_err_bssCSD_opt.mat')
load('/export/home/mike/matlab/batch/ac252_err_bssCSD_fpo93/ac252_err_bssCSD_fpo93_opt.mat')

%{
opts = { '/export/home/mike/matlab/batch/ac252_vp3_nonCSD_mono/ac252_vp3_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_cpt_nonCSD_mono/ac252_cpt_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_ern_nonCSD_mono/ac252_ern_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_aod_nonCSD_mono/ac252_aod_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_ans_nonCSD_mono/ac252_ans_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_stp_nonCSD_mono/ac252_stp_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_gng_nonCSD_mono/ac252_gng_nonCSD_mono_opt.mat'};
%}

%{
opts = {'/export/home/mike/matlab/batch/ac252_vp3_noasr_mono/ac252_vp3ac252_vp3_noasr_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_vp3_asr_mono/ac252_vp3ac252_vp3_asr_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_vp3_asr_noburst/ac252_vp3ac252_vp3_asr_noburst_opt.mat'};
%}

%{
opts = {'/export/home/mike/matlab/batch/nki_vp3_noclean_noeog_simple/nki_vp3_cnth1s_noclean_noeog_simple_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_vp3_clean_noeog_simple/nki_vp3_cnth1s_clean_noeog_simple_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_vp3_clean_eog_simple/nki_vp3_cnth1s_clean_eog_simple_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_vp3_clean_crlseog_simple/nki_vp3_cnth1s_clean_crlseog_simple_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_vp3_crlseog_both/nki_vp3_cnth1s_crlseog_both_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_vp3_bsseog_both/nki_vp3_cnth1s_bsseog_both_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_vp3_bsseog_newboth/nki_vp3_cnth1s_bsseog_newboth_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_vp3_bsseognew_both/nki_vp3_cnth1s_bsseognew_both_opt.mat'
    };
%}

%{
opts = {'/export/home/mike/matlab/batch/nki_err_nonCSD/nki_err_cnth1snonCSD_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_err_bsseog_newboth/nki_err_cnth1s_bsseog_newboth_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_err_bssnew_nonCSD/nki_err_cnth1s_bssnew_nonCSD_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_err_CSD/nki_err_cnth1sCSD_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_err_bsseog_newboth_CSD/nki_err_cnth1s_bsseog_newboth_CSD_opt.mat', ...
    '/export/home/mike/matlab/batch/nki_err_bssnew_CSD/nki_err_cnth1s_bssnew_CSD_opt.mat'};
%}

if exist('opts', 'var')
    opt = coh_multiopt( opts );
end

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
%demogsfile='/export/home/mike/matlab/batch/err_csub_FPO90/err_csubs_csub_FPO90_demog.mat'; %***
%demogsfile='/export/home/mike/matlab/batch/asubs_err_fpo90/err_asubs_FPO90_demog.mat';
%demogsfile='/export/home/mike/matlab/batch/asubs_err_fpo90/err_acsubscombined_FPO90_demog302_6col.mat'; % ****
%demogsfile='/export/home/mike/matlab/batch/asubs_err_fpo90/err_acsubscombined_FPO90_demog160_6col.mat';
%demogsfile='/export/home/mike/matlab/csv/ashwini_38.mat';
%demogsfile='/export/home/mike/matlab/csv/kam_err_ph4_ingroup_noiowa.mat';
%demogsfile='/export/home/mike/matlab/batch/fmri_err_FPO90/fmri_err_FPO90_demogs.mat';
%demogsfile='/export/home/mike/matlab/csv/nki.mat';
%demogsfile=[];
%demogsfile='/export/home/mike/matlab/csv/mf90.mat';
demogsfile='/export/home/mike/matlab/csv/mf90_42.mat';

% set the factor for time-downsampling
time_ds=4;

% set the path for behavioral data if applicable
%behmat_path='/active_projects/mike/fmri_phase4_ern_beh';
%behmat_path='/active_projects/mike/phase4_ern_beh';
behmat_path='/active_projects/mike/ern_csub_beh';
%behmat_path=[];

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
datatypes={'n_trials','erp','wave_totpow','wave_evknorm','coh'};
%datatypes={'n_trials','erp','wave_totpow','wave_evknorm'};
%datatypes={'n_trials','erp'};

%[mat_list, imp, n_trials_all, behdata, erpdata, ~, wave_evknormdata, cohdata, phidata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path, datatypes);
% LONG ERP
%[mat_list, imp, n_trials_all, behdata, erpdata, wave_totpowdata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path, datatypes);
% REGULAR
%[mat_list, imp, n_trials_all, behdata, erpdata, ~, wave_evknormdata, cohdata, phidata] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path, datatypes);
% KAM / ASHWINI
%[mat_list, imp, n_trials_all, ~, erpdata, wave_totpowdata, wave_evknormdata, ~, ~, ~] = ...
%    coh_import(opt, demogsfile, time_ds, behmat_path, datatypes);

% multi-opt
if length(opt) ~= 1
    [mat_list, imp, n_trials_all, ~, erpdata, ~, ~, ~, ~, ~] = ...
        coh_multimport(opt, demogsfile, time_ds, behmat_path, datatypes);
else
    [mat_list, imp, n_trials_all, behdata, erpdata, wave_totpowdata, wave_evknormdata, cohdata, ~, ~] = ...
        coh_import(opt, demogsfile, time_ds, behmat_path, datatypes);
end

%%

% scales
[scl,chan_locs] = coh_scales(opt, imp, -187.5);

% rejection & demographics
trials_necessary=0; %15 for ERP/power, 40 for ITC/ISPC
%custom_rej=[7 10 14 23 28 41 52 54 57];
%custom_rej=[12,28,38,39,41,43,61,63,65,68,82,86,87,89,97,99,101,108,115,117,119,124,126,129,133,134,135,146,157,159,160,161];
%custom_rej=[1,5,6,8,14,22,23,26,27,33,40,41,45,47,54,55,58,60,64,69,78,82,93,102,103,110,115,118,127,147,158,167,168,169,175,182,189,194]; %csubs
%[230,242,540,545,1044,1345,1457,1464,1468,1483,1492,1825,1826,1839,2030,2034]; %kam
%custom_rej = [5,6,13,17,19,30,36,43,44,48,59,60,61,62,82,84,93,97,110,112,114,116,117,118,119,121,153,158,166,171,173]; %asubs
%custom_rej = [3,6,12,15,21,53,55,57,58,60,87,90,94,98,101,109,111,135,178]; %acsubs
custom_rej=[];
%age_range=[21 45];
age_range=[];
[scl, s_inds_g, s_demogs] = coh_demog(mat_list, imp, scl, n_trials_all, demogsfile, ...
    trials_necessary, age_range, custom_rej);

% plotting parameters
cond_diff = {1, 2};
[scl, pp] = coh_plotparams(imp, scl, cond_diff);

%%

% filter / baseline / pick peaks in ERP
%cutoff = [2 16];
cutoff = 16;
[erpdata,peakmat] = coh_pickERPpeaks(imp, scl, s_inds_g, erpdata, [true true], cutoff);

%%

clear datatypes demogsfile time_ds behmat_path lp_cutoff custom_rej ...
    trials_necessary age_range cond_diff