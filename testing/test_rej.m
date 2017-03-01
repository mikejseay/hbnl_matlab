%% hackish way to load data (if so, skip down)

load('/active_projects/mike/test_trialdata/ern_s4.mat')
artf_chan_vec = 1:61;

%% load cleaned data instead

start_dir = '/processed_data/matlab_common/ern/clean/';
files = dir(start_dir);
files = {files.name}';
files(1:2) = [];

load('/export/home/mike/matlab/batch/nki_ern_nonCSD/nki_ern_cnth1snonCSD_opt.mat')

%% hack to fix the fact that the memoized cleaned/CSD data does not save the name of the h1 file

filecell = list2cell(opt.infile_list);
n_files = length(filecell);
h1_uids = cellfun( @(x)x(63:73), filecell, 'uni', 0);
mat_uids = cellfun( @(x)x(7:17), files, 'uni', 0);

[~, h1_sort] = sort(h1_uids);
[~, mat_sort] = sort(mat_uids);

h1_files = filecell(h1_sort);
mat_files = files(mat_sort);

%% choose a file to load

opt.simple_thresh = 150;
opt.eeglab_thresh = 75;

prop_kept_simple_all = zeros(60,1);
prop_kept_eeglab_all = zeros(60,1);
prop_kept_both_all = zeros(60,1);
trial_mat_simple_all = cell(60,1);
trial_mat_eeglab_all = cell(60,1);
trial_mat_both_all = cell(60,1);

for file=[8:60]

h1_struct = read_hdf1_dataV7(h1_files{file});
data_file = files{file};
load(fullfile(start_dir, data_file))
[~, n_samps, dataR, trial_type] = hdf_reshape(dataR, h1_struct, opt, trial_type);


n_trials = size(dataR, 2);
n_cases = max(trial_type);
% if no custom trial indexings were done above
if ~isfield(opt, 'trial_mat')

    % make trial index based on type codes
    trial_mat = false(n_trials, n_cases);
    for m = 1:n_cases
        trial_mat(:, m) = trial_type == m;
    end
    if isfield(opt, 'case_vec') %in other words, if there is a
        % case vector specified, take only those, otherwise take all
        % cases
        trial_mat = trial_mat(:, opt.case_vec);
    end
    zero_trials = find(sum(trial_mat) == 0);
    if ~isempty(zero_trials)
        trial_mat(:, zero_trials) = [];
    end
else
    trial_mat = opt.trial_mat;
end
n_cases = size(trial_mat, 2);

% trial rejection

% first specify the channel / timepoints range to examine
if isfield(opt, 'artf_chan_vec')
    artf_chan_vec = opt.artf_chan_vec; %use specified channels if specified
else
    artf_chan_vec = 1:n_chans; %or all of them if not
end

if isfield(opt, 'artf_samp_vec')
    artf_samp_vec = opt.artf_samp_vec; %use certain samples if specified
elseif isfield(opt, 'artf_ms_vec')
    artf_samp_vec = round(opt.artf_ms_vec/1000*...
        opt.rate);
    artf_samp_vec=artf_samp_vec(1):artf_samp_vec(2); %or a ms range, better
end
if isfield(opt, 'prestim_ms')
    artf_samp_vec = artf_samp_vec + round(opt.prestim_ms/1000 * opt.rate) + 1;
elseif isfield(opt, 'epoch_lims')
    artf_samp_vec = artf_samp_vec + round(abs(opt.epoch_lims(1))/1000 * opt.rate) + 1;
else
    artf_samp_vec = 49:166;
end

h1_struct = struct();
h1_struct.experiment_struct.rate = opt.rate; %hack to make eeglab rejection work below

disp(['done loading',num2str(file)]);
trial_mat_orig = trial_mat;

%

% simple

trial_mat_simple = check_artifact(dataR(artf_samp_vec, :, :), trial_mat_orig, ...
        artf_chan_vec, opt.simple_thresh);
prop_kept_simple=sum(sum(trial_mat_simple))/sum(sum(trial_mat_orig));
prop_kept_simple_all(file) = prop_kept_simple;
trial_mat_simple_all{file} = trial_mat_simple;

%

% eeglab

for sdthresh=3.5:.25:5
    trial_mat_eeglab = check_artifact_eeglab(data_file, dataR(artf_samp_vec, :, :), ...
        h1_struct, trial_mat_orig, sdthresh, opt.eeglab_thresh, artf_chan_vec);
    prop_kept_eeglab=sum(sum(trial_mat_eeglab))/sum(sum(trial_mat_orig));
    if prop_kept_eeglab >= .7
        fprintf(1, 'SD threshold used was %1.2f\n', sdthresh);
        break
    end
end
prop_kept_eeglab_all(file) = prop_kept_eeglab;
trial_mat_eeglab_all{file} = trial_mat_eeglab;

%

% both

trial_mat_both = repmat(trial_mat, [1 1 2]);
trial_mat_both(:,:,1) = trial_mat_simple;
trial_mat_both(:,:,2) = trial_mat_eeglab;
trial_mat_both2 = squeeze(all(trial_mat_both, 3));
prop_kept_both = sum(sum(trial_mat_both2))/sum(sum(trial_mat_orig));
%if prop_kept_both <= .65
%    disp('simple thresh was overly aggressive, taking only eeglab rejects')
%    %trial_mat = squeeze(any(trial_mat_both, 3)); %experimental
%    trial_mat_both2 = trial_mat_eeglab;
%end

sum(trial_mat_both(:,:,1))
sum(trial_mat_both(:,:,2))

prop_kept_both_all(file) = prop_kept_both;
trial_mat_both_all{file} = trial_mat_both;

fprintf('file %d - simple: %0.2f, eeglab: %0.2f, both: %0.2f\n', ...
    file, prop_kept_simple, prop_kept_eeglab, prop_kept_both);

end

%% filtered/baselined threshold

opt.thresh_bl = [-100 0];

if isfield(opt, 'thresh_bl') && isfield(opt, 'artf_ms_vec')
    ms_vec = opt.artf_ms_vec(1):1000/opt.rate:opt.artf_ms_vec(2);
    ms_vec(end) = [];
    bl_pts = dsearchn(ms_vec', [-100, 0]');
    bl_vec = bl_pts(1):bl_pts(end);
else
    bl_vec = 103:129;
end
trial_mat = check_artifact_baselined(dataR(artf_samp_vec, :, :), trial_mat_orig, ...
    artf_chan_vec, opt.thresh, bl_vec);