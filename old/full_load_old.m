ern_cleannewconds_ql %specify the path, demographics file, and chan_locs
%coh_import %import selected data types
coh_import_older
coh_scales_old %scale dimensions
coh_reject_old %reject subjects based on trials
coh_demog_old %match demographics and index subjects
coh_plotparams_old %specify plotting parameters
coh_pickERPpeaks_old %filter/baseline ERPs and pick peaks in a time window

%load('/export/home/mike/matlab/origin/fmri/fmri_ern_beh.mat') %load
%behavioral data
