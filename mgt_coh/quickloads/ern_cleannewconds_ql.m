%% mat path (=opt.outpath), demographics file (must be supplied)

%matpath='/active_projects/mike/fmri_phase4_ern_newconds/';
%matpath='/active_projects/mike/fmri_phase4_err_noclean/';
%matpath='/active_projects/mike/pcheck_condmean/';
%matpath='/active_projects/mike/fmri_phase4_err_noclean36/';
%matpath='/active_projects/mike/fmri_phase4_err_clean36/';
%matpath='/active_projects/mike/fmri_phase4_err_cleanCSD/';
%matpath='/active_projects/mike/fmri_phase4_err_nocleanfilt/';
%matpath='/active_projects/mike/fmri_phase4_err_cleanCSD3seeds/';
%matpath='/active_projects/mike/fmri_phase4_err_clean_hpfilttest/';
%matpath='/active_projects/mike/fmri_phase4_err_clean_hpfilttest2/';
%matpath='/active_projects/mike/fmri_phase4_err_cleanCSD84/';
%matpath='/active_projects/mike/fmri_phase4_ern_cleanCSD_newconds2/';
%matpath='/active_projects/mike/fmri_phase4_ern_cleanCSD4seeds/';
%matpath='/active_projects/mike/fmri_phase4_err_cleanCSD2seeds/';
%matpath='/active_projects/mike/fmri_phase4_ern_beh/';
%matpath='/active_projects/mike/fmri_phase4_err_cleanCSD2seeds_fpzoz/';
%matpath='/active_projects/mike/fmri_phase4_err_cleanCSD4seeds/';
%matpath='/active_projects/mike/fmri_phase4_err_cleanCSD_eprej_4seeds/';
matpath='C:\Users\Mike\Desktop\localdata\fmri_err_FPO92';
%matpath='/active_projects/mike/fmri_err_FPO92/';

demogsfile='A:/matlab/origin/fmri/fmri_demogs_simple.mat';

if true
    mat_list=coh_updatemats(matpath,demogsfile);
else
    mat_list=coh_updatemats(matpath);
end

%% opt

%load('/export/home/mike/matlab/batch/fmri_ern_wavelet_newconds.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_noclean.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_clean.mat')
%load('/export/home/mike/matlab/batch/ern_fmri_pchecknew.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_noclean36.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_clean36.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_cleanCSD.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_nocleanfilt.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_cleanCSD3seeds.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_clean_hpfilttest.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_clean_hpfilttest2.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_cleanCSD84.mat')
%load('/export/home/mike/matlab/batch/fmri_ern_wavelet_cleanCSD_newconds2.mat')
%load('/export/home/mike/matlab/batch/fmri_ern_wavelet_cleanCSD4seeds.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_cleanCSD2seeds.mat')
%load('/export/home/mike/matlab/batch/fmri_ern_beh.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_cleanCSD4seeds.mat')
%load('/export/home/mike/matlab/batch/fmri_err_wavelet_cleanCSD_eprej_4seeds.mat')
%load('A:/matlab/batch/err_wavelet_cleanCSD_eprej_4seeds.mat');
load('A:/matlab/batch/err_FPO92/fmri_err_cnth1s_FPO92_opt.mat')

if exist('param_struct','var')
    opt=fix_opt(param_struct);
end


%% chan_locs

load('A:/matlab/origin/coords/61chans_ns.mat')

clear matpath demogsfile