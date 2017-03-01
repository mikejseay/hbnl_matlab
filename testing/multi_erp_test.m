% determine what is the set intersection of the subjects available for a
% set of folders containing .mat files, starting from their opts
% vp3, cpt, ern, aod, ans, stp, gng

demogsfile='/export/home/mike/matlab/batch/asubs_err_fpo90/err_acsubscombined_FPO90_demog302_6col.mat';

opts = { '/export/home/mike/matlab/batch/ac252_vp3_nonCSD_mono/ac252_vp3_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_cpt_nonCSD_mono/ac252_cpt_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_ern_nonCSD_mono/ac252_ern_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_aod_nonCSD_mono/ac252_aod_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_ans_nonCSD_mono/ac252_ans_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_stp_nonCSD_mono/ac252_stp_nonCSD_mono_opt.mat', ...
    '/export/home/mike/matlab/batch/ac252_gng_nonCSD_mono/ac252_gng_nonCSD_mono_opt.mat'};
n_opts = length(opts);

matlist_all = cell(n_opts, 2);

for o = 1:n_opts
    load( opts{o} )
    [matlist_all{o, 2}, matlist_all{o, 1}] = coh_updatemats(opt.outpath, demogsfile);
end

common_ids = mintersect(matlist_all(:,1));

%% 

% create a modified demogsfile that only uses the common ids

load(demogsfile)
demogs_commonids = cell2table(common_ids, 'VariableNames', {'UniqueFileString'});
demogs_table = join(demogs_commonids, demogs_table, 'Keys', [1 1]);