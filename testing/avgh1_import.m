%% constants
chan_vec = [1:31,33:62];
sen = 204.8./6.1035156;

%% params

paths = {'/export/home/mike/avg_hp003_t75/', '/export/home/mike/avg_hp05_t150/'};
demogsfile = '/export/home/mike/matlab/csv/nki.mat';

%% load

erpdata_h1 = zeros(416, 61, 6, 60);
conds = 1:3;
for p=1:length(paths)
    path = paths{p};
    avgh1_list = coh_updatemats(path,demogsfile,'.h1');

    for f=1:length(avgh1_list)

        h1_path = fullfile(avgh1_list{f});
        h1_struct = read_hdf1_dataV7(h1_path);
        data = permute(h1_struct.data_struct.hdf1_avg_data, [3 2 1]);
        erpdata_h1(:, :, conds, f) = data(:, chan_vec, :).*sen;
        clear data

    end
    conds = conds + 3;
end

%% mocking up scaling and plotting assist structures

scl_h1 = struct();
scl_h1.t_ms_erp = -h1_struct.experiment_struct.pre_stim_time_ms:(1000/256): ...
    h1_struct.experiment_struct.post_stim_time_ms;
scl_h1.t_ms_erp(end) = [];
scl_h1.erpmaxtimepts = length(scl_h1.t_ms_erp);
scl_h1.t_xtick_erp = dsearchn(scl_h1.t_ms_erp'/1000, scl.t_xtick_ms');
scl_h1.t_zero_erp = dsearchn(scl_h1.t_ms_erp'/1000, 0);
[~, scl_h1.t_start_ind_erp] = min(abs(scl_h1.t_ms_erp - scl.t_ms_erp(scl.t_start_ind)));
[~, scl_h1.t_end_ind_erp] = min(abs(scl_h1.t_ms_erp - scl.t_ms_erp(imp.erpmaxtimepts)));

scl_h1.cond_label([1 4]) = {'avgh1-003-t75', 'avgh1-05-t150'};

pp_h1.plotn_cond = [1 4];
pp_h1.cond_diff = {[1],[2]};
pp_h1.sp_d = numSubplots(length(pp_h1.plotn_cond));

%% clean up

clear path demogsfile chan_vec sen h1_path h1_struct data conds