%% load the data

load('higherror_workspace.mat');
time_vec = double(time_vec);

%% original scheme which is just a strict downsampling

f_ds_vec = (f>=0.5 & f<=45)';
f_vec = f';
f_vec_ds = f(f_ds_vec)';

rate = h1_struct.experiment_struct.rate;
time_ds_factor = 4;
time_vec = (-h1_struct.experiment_struct.pre_stim_time_ms:(1000/rate): ...
    h1_struct.experiment_struct.post_stim_time_ms)';
time_cut_vec = (time_vec>=-200 & time_vec<=900);
time_downsamp_vec = false(size(time_vec));
time_downsamp_vec(1:time_ds_factor:end) = true;
time_ds_vec = time_cut_vec & time_downsamp_vec;
time_vec_ds = time_vec(time_ds_vec);

opt = v2struct(add_baseline, calc_type, case_name, channel_sort, do_baseline, ...
    elec_array, exp_name, file_id, filenm, file_run, ...
    file_session, i_file, ln_calc, n_chans, out_type, out_type_name, ...
    output_text, st_type, st_type_name, ...
    f_vec, f_vec_ds, rate, time_ds_factor, time_vec, time_vec_ds);

S_all_ds = S_all(f_ds_vec, time_ds_vec, :);
% clear S_all
data = single(S_all_ds);
% clear S_all_ds

folderparts = strsplit(filenm,'/');
nameparts = strsplit(folderparts{end},'_');
session_part = nameparts{3};
id = nameparts{4};
session = session_part(1);
param_str = folderparts{end-2};

parent_dir = '/processed_data/ero-mats';
outname = [strjoin({id, session, exp_name, case_name, st_type_name}, '_'),'.mat'];
outdir = fullfile(parent_dir, param_str, num2str(n_chans), exp_name);
outpath = fullfile(outdir, outname);

%% full data (no downsampling)
% recall the crucial fields are f_vec_ds and time_vec_ds
% and the data obj should be called data

data = single(S_all);
time_vec_ds = time_vec;
f_vec_ds = f_vec;
opt = v2struct(add_baseline, calc_type, case_name, channel_sort, do_baseline, ...
    elec_array, exp_name, file_id, filenm, file_run, ...
    file_session, i_file, ln_calc, n_chans, out_type, out_type_name, ...
    output_text, st_type, st_type_name, ...
    f_vec, f_vec_ds, rate, time_ds_factor, time_vec, time_vec_ds);

outdir = '/export/home/mike/matlab/test_eromats';
outname = [strjoin({id, session, exp_name, case_name, st_type_name}, '_')];
outname = [outname, '_full'];
outpath = fullfile(outdir, [outname, '.mat']);
save(outpath, 'data', 'opt', '-v7.3');

time_win = [200 400];
freq_win = [4 8];

[~, t_start] = min(abs(time_vec - time_win(1)));
[~, t_end] = min(abs(time_vec - time_win(2)));
[~, f_start] = min(abs(f_vec - freq_win(1)));
[~, f_end] = min(abs(f_vec - freq_win(2)));

tf_means = squeeze(mean(mean(S_all(f_start:f_end, t_start:t_end, :), 1), 2));

%% downsampled the way i initially designed it

f_ds_vec = (f>=0.5 & f<=45)';
f_vec = f';
f_vec_ds = f(f_ds_vec)';

rate = h1_struct.experiment_struct.rate;
time_ds_factor = 1;
time_vec = (-h1_struct.experiment_struct.pre_stim_time_ms:(1000/rate): ...
    h1_struct.experiment_struct.post_stim_time_ms)';
time_cut_vec = (time_vec>=-200 & time_vec<=900);
time_downsamp_vec = false(size(time_vec));
time_downsamp_vec(1:time_ds_factor:end) = true;
time_ds_vec = time_cut_vec & time_downsamp_vec;
time_vec_ds = time_vec(time_ds_vec);

opt = v2struct(add_baseline, calc_type, case_name, channel_sort, do_baseline, ...
    elec_array, exp_name, file_id, filenm, file_run, ...
    file_session, i_file, ln_calc, n_chans, out_type, out_type_name, ...
    output_text, st_type, st_type_name, ...
    f_vec, f_vec_ds, rate, time_ds_factor, time_vec, time_vec_ds);

S_all_ds = S_all(f_ds_vec, time_ds_vec, :);
% clear S_all
data = single(S_all_ds);
% clear S_all_ds

outdir = '/export/home/mike/matlab/test_eromats';
outname = [strjoin({id, session, exp_name, case_name, st_type_name}, '_')];
outname = [outname, '_ds', num2str(time_ds_factor)];
outpath = fullfile(outdir, [outname, '.mat']);
save(outpath, 'data', 'opt', '-v7.3');

time_win = [200 400];
freq_win = [4 8];

[~, t_start] = min(abs(time_vec_ds - time_win(1)));
[~, t_end] = min(abs(time_vec_ds - time_win(2)));
[~, f_start] = min(abs(f_vec_ds - freq_win(1)));
[~, f_end] = min(abs(f_vec_ds - freq_win(2)));

tf_means = squeeze(mean(mean(data(f_start:f_end, t_start:t_end, :), 1), 2));

%% downsample with average

f_ds_vec = (f>=0.5 & f<=45)';
f_vec = f';
f_vec_ds = f(f_ds_vec)';

rate = h1_struct.experiment_struct.rate;
time_ds_factor = 4;
time_vec = (-h1_struct.experiment_struct.pre_stim_time_ms:(1000/rate): ...
    h1_struct.experiment_struct.post_stim_time_ms)';
time_cut_vec = (time_vec>=-200 & time_vec<=900);
time_downsamp_vec = false(size(time_vec));
time_downsamp_vec(1:time_ds_factor:end) = true;
time_ds_vec = time_cut_vec & time_downsamp_vec;
time_vec_ds = time_vec(time_ds_vec);

S_all_fds = S_all(f_ds_vec, :, :);
S_all_ds = zeros(size(S_all_fds, 1), length(time_vec_ds), size(S_all_fds, 3));
time_inds = find(time_ds_vec);
time_vec_ds_avg = zeros(length(time_inds), 1);
for t=1:length(time_inds)
    S_all_ds(:, t, :) = squeeze(mean(S_all_fds(:, time_inds(t):time_inds(t)+3, :), 2));
    time_vec_ds_avg(t) = mean(time_vec(time_inds(t):time_inds(t)+3));
end
% clear S_all
data = single(S_all_ds);
% clear S_all_ds

time_vec_ds = time_vec_ds_avg;
opt = v2struct(add_baseline, calc_type, case_name, channel_sort, do_baseline, ...
    elec_array, exp_name, file_id, filenm, file_run, ...
    file_session, i_file, ln_calc, n_chans, out_type, out_type_name, ...
    output_text, st_type, st_type_name, ...
    f_vec, f_vec_ds, rate, time_ds_factor, time_vec, time_vec_ds);

outdir = '/export/home/mike/matlab/test_eromats';
outname = [strjoin({id, session, exp_name, case_name, st_type_name}, '_')];
outname = [outname, '_ds4avg_newind'];
outpath = fullfile(outdir, [outname, '.mat']);
save(outpath, 'data', 'opt', '-v7.3');

time_win = [200 400];
freq_win = [4 8];

[~, t_start] = min(abs(time_vec_ds - time_win(1)));
[~, t_end] = min(abs(time_vec_ds - time_win(2)));
[~, f_start] = min(abs(f_vec_ds - freq_win(1)));
[~, f_end] = min(abs(f_vec_ds - freq_win(2)));

tf_means = squeeze(mean(mean(data(f_start:f_end, t_start:t_end, :), 1), 2));

%% using interp2 

time_vec = double(time_vec);

time_start = time_vec(1);
time_space = 12.5;
time_end = 900;

freq_start = 0.5;
freq_space = 0.5;
freq_end = 45;

time_vec_interp = time_start:time_space:time_end;
freq_vec_interp = freq_start:freq_space:freq_end;
[tq, fq] = meshgrid(time_vec_interp, freq_vec_interp);
n_times_interp = size(tq, 1);
n_freqs_interp = size(tq, 2);
S_all_it2 = zeros(n_times_interp, n_freqs_interp, n_chans);
S_all_it2_lin = zeros(n_times_interp, n_freqs_interp, n_chans);
S_all_it2_cub = zeros(n_times_interp, n_freqs_interp, n_chans);
S_all_it2_spl = zeros(n_times_interp, n_freqs_interp, n_chans);

% (linear, default)

for chan=1:n_chans
    S_all_it2_lin(:, :, chan) = interp2(time_vec, f_vec, squeeze(S_all(:, :, chan)), tq, fq);
end

% (cubic)

for chan=1:n_chans
    S_all_it2_cub(:, :, chan) = interp2(time_vec, f_vec, squeeze(S_all(:, :, chan)), tq, fq, 'cubic');
end

% (spline)

for chan=1:n_chans
    S_all_it2_spl(:, :, chan) = interp2(time_vec, f_vec, squeeze(S_all(:, :, chan)), tq, fq, 'spline');
end

%% linear interp save

time_vec_ds = time_vec_interp';
f_vec_ds = freq_vec_interp';

data = single(S_all_it2_lin);

opt = v2struct(add_baseline, calc_type, case_name, channel_sort, do_baseline, ...
    elec_array, exp_name, file_id, filenm, file_run, ...
    file_session, i_file, ln_calc, n_chans, out_type, out_type_name, ...
    output_text, st_type, st_type_name, ...
    f_vec, f_vec_ds, rate, time_ds_factor, time_vec, time_vec_ds);

outdir = '/export/home/mike/matlab/test_eromats';
outname = [strjoin({id, session, exp_name, case_name, st_type_name}, '_')];
outname = [outname, '_interplin'];
outpath = fullfile(outdir, [outname, '.mat']);
save(outpath, 'data', 'opt', '-v7.3');

time_win = [200 400];
freq_win = [4 8];

[~, t_start] = min(abs(time_vec_ds - time_win(1)));
[~, t_end] = min(abs(time_vec_ds - time_win(2)));
[~, f_start] = min(abs(f_vec_ds - freq_win(1)));
[~, f_end] = min(abs(f_vec_ds - freq_win(2)));

tf_means = squeeze(mean(mean(data(f_start:f_end, t_start:t_end, :), 1), 2));

%% plot all three

chan=25;

figure;
subplot(131); contourf(S_all_it2_lin(:, :, chan));
subplot(132); contourf(S_all_it2_cub(:, :, chan));
subplot(133); contourf(S_all_it2_spl(:, :, chan));

%% using interp1 

% interp2 seems strictly better?

%% using imresize

% this requires the image processing toolbox and should be avoided
