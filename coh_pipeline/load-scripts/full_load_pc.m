% new full load

% load the options file (as opt here)
% (or param_struct if legacy)
%load('A:/matlab/batch/err_FPO92/fmri_err_cnth1s_FPO92_opt.mat')
load('C:\Users\Mike\Documents\MATLAB\batch\fmri_err_wavelet_cleanCSD4seeds.mat')

if exist('param_struct','var')
    outpath = 'C:\Users\Mike\Desktop\localdata\fmri_phase4_err_cleanCSD4seeds';
    coordsfile = 'C:\Users\Mike\Documents\MATLAB\origin\coords\61chans_ns.mat';
    opt=fix_opt(param_struct,outpath,coordsfile);
    clear param_struct outpath coordsfile
end

% set the demographics file if applicable
demogsfile='A:/matlab/origin/fmri/fmri_demogs_simple.mat';

% set the factor for time-downsampling
time_ds=1;

% import
%[mat_list, imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata] = ...
%    coh_import(opt, demogsfile);
[mat_list, imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata] = ...
    coh_import_legacy(opt, demogsfile, time_ds);

% scales
[scl,chan_locs] = coh_scales(opt,imp);

% rejection & demographics
[scl, s_inds_g, s_demogs] = coh_demog(mat_list, imp, scl, n_trials_all, demogsfile);

% plotting parameters
[scl,pp] = coh_plotparams(imp,scl);

% filter / baseline / pick peaks in ERP
lp_cutoff = 24;
[erpdata,peakmat] = coh_pickERPpeaks(imp, scl, s_inds_g, erpdata, lp_cutoff);

clear demogsfile time_ds lp_cutoff