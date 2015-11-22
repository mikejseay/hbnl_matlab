matpath='/active_projects/mike/fmri_phase4_err_cleanCSD/';
demogsfile='/export/home/mike/matlab/origin/fmri/fmri_demogs_simple.mat';
mat_list=coh_updatemats(matpath,demogsfile);
load('/export/home/mike/matlab/batch/fmri_err_wavelet_cleanCSD.mat')
load('/export/home/mike/matlab/origin/coords/61chans_ns.mat')
clear matpath demogsfile