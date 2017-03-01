% import all CNTs in a folder together

% get expstruct that defines experiment properties
build_expstruct;

% define folder to target
target_dir = '/raw_data/neuroscan/suny/ns620/49432071';

% define set of experiments to extract
exps = {'vp3'};
%exps = {}; % empty defaults to all

% set whether to merge them
merge_bool = true;

% define epoch limits
ep_lims = [-.4 .8];
%ep_lims = [-.65 1];

% some pre-defined information
nonevent_exps = {'eec', 'eeo'};
remove_chans_str = {'BLANK'};
hp_cutoff = 1;
resample_rate = 250;
nonhead_chans = [32 63];
nonhead_chans_str = {'X', 'Y'};
%head_chans = [1:31, 33:62];
head_chans_str = {'FP1','FP2','F7','F8','AF1','AF2','FZ','F4','F3','FC6', ...
    'FC5','FC2','FC1','T8','T7','CZ','C3','C4','CP5','CP6','CP1','CP2', ...
    'P3','P4','PZ','P8','P7','PO2','PO1','O2','O1','AF7','AF8','F5','F6', ...
    'FT7','FT8','FPZ','FC4','FC3','C6','C5','F2','F1','TP8','TP7','AFZ', ...
    'CP3','CP4','P5','P6','C1','C2','PO7','PO8','FCZ','POZ','OZ','P2','P1','CPZ'};
head_chans = 1:61;
%head_eye_chans = [1 2 5 6 33 34 39 48];
head_eye_chans = [1 2 5 6 32 33 38 47];
head_eye_chans_str = {'FP1','FP2','AF1','AF2','AF7','AF8','FPZ','AFZ'};
head_noneye_chans = setdiff(head_chans, head_eye_chans);
head_noneye_chans_str = setdiff(head_chans_str, head_eye_chans_str);
uvthresh = 100;
sdthresh = 4;
%max_rejperecent = 15;

% index its contents, pull out names only
contents = dir(target_dir);
file_names = { contents.name }';

% find .cnt files, extract their experiment names
% exclude if name contains "orig"
% sort in recording order (always same)
cnt_inds = find( ~cellfun('isempty',regexp(file_names, 'cnt')) & ...
    cellfun('isempty',regexp(file_names, 'orig')) );
dat_inds = find( ~cellfun('isempty', regexp(file_names, 'dat')) );
cnt_names = file_names( cnt_inds );
dat_names = file_names( dat_inds );
experiment_names_cnt = cellfun( @(x)x(1:3), cnt_names, 'uni', 0);
experiment_names_dat = cellfun( @(x)x(1:3), dat_names, 'uni', 0);

% get indices of experiment to target
if isempty(exps)
    exp_ind_set = [1:length(expstruct)]';
else
    [~, exp_ind_set] = intersect({expstruct.name}, exps);
    exp_ind_set = sort(exp_ind_set);
end

%%

% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% import each experiment, renaming trials
for exp_ind = exp_ind_set
    
    % if that experiment is missing, skip it
    if ~ismember(expstruct(exp_ind).name, experiment_names_cnt)
        continue
    end
    
    cfile = find( strcmpi(experiment_names_cnt, expstruct(exp_ind).name) );
    dfile = find( strcmpi(experiment_names_dat, expstruct(exp_ind).name) );
    
    target_cntfile = cnt_names{ cfile };
    target_datfile = dat_names{ dfile };
    target_cntpath = fullfile( target_dir, target_cntfile );
    target_datpath = fullfile( target_dir, target_datfile );
    
    % load the CNT
    EEG = pop_loadcnt(target_cntpath, 'dataformat', 'auto', 'memmapfile', '');
    
    % import behavioral info by attempting to access the same-named
    % *.dat file in the same folder
    %[~, dat_name] = fileparts(target_file);
    %dat_file = [dat_name, '.dat'];
    %target_datpath = fullfile(target_dir, dat_file);
    try
        EEG = pop_loaddat(EEG, target_datpath, expstruct(exp_ind).iti);
    catch
        fprintf('No dat present for %s\n',expstruct(exp_ind).name);
    end
    
    % rename the trials
    EEG = rename_EEG(EEG, expstruct(exp_ind).ttl_maps);
    
    % if one of the no-event conditions, add generic markers every epoch
    % length
    if ismember(experiment_names_cnt{cfile}, nonevent_exps)

        [EEG.event, EEG.urevent] = build_evstruct( ...
            length(EEG.times), EEG.srate, sum(abs(ep_lims)), ...
            [experiment_names_cnt{cfile},'_seg'] );
        
    end
    
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, ...
        cfile - 1, 'setname', target_cntfile(1:end-4), 'gui', 'off');

end

% merge
if merge_bool && length(ALLEEG) > 1
    EEG = eeg_checkset( EEG );
    EEG = pop_mergeset( ALLEEG, 1:length(ALLEEG), 0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, length(cnt_names), 'gui', 'off'); 
end

% delete remainder (?)
%[ALLEEG EEG CURRENTSET] = pop_copyset(ALLEEG, length(ALLEEG), 1 );
%ALLEEG = pop_delset( ALLEEG, 2:length(ALLEEG) );
%ALLEEG(2:end) = [];

%%

% delete blank channel
EEG = eeg_checkset( EEG );
EEG = pop_select( EEG, 'nochannel', remove_chans_str);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off'); 

% import the channel locations
load('/export/home/mike/matlab/origin/coords/63chans_ns_eog.mat')
EEG.chanlocs=chan_locs;

% delete non-head channels (?)
EEG = eeg_checkset( EEG );
EEG = pop_select( EEG, 'nochannel', nonhead_chans_str);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off'); 

% use clean_channels (excluded for now)
%{
% with only channel rejection (?), no segment rejection
% also turned off channel correlation and line noise method because it's slow
EEG = eeg_checkset( EEG );
n_chans = size(EEG.data,1);
EEG = cleanset_eeglab(EEG, [0.25 0.75], 'off', 'off', 'off');

% record the names of and interpolate any channels removed by clean_channels
n_interpchans = n_chans - size(EEG.data,1);
if n_interpchans > 0
    EEG = pop_interp(EEG, chan_locs, 'spherical');
end
%}

% hi-pass filter at 1 Hz (only if didn't use clean_channels above)
EEG = eeg_checkset( EEG );
EEG = pop_eegfiltnew(EEG, hp_cutoff, []);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off'); 

% resample to 250 Hz
EEG = eeg_checkset( EEG );
EEG = pop_resample( EEG, resample_rate);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off'); 

% re-reference head electrodes to average of head electrodes
EEG = eeg_checkset( EEG );
%EEG = pop_reref( EEG, [],'exclude', nonhead_chans);
EEG = pop_reref( EEG, []);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off'); 

% epoch based on defined limits (e.g. -.4 to .8)
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  }, ep_lims, 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off'); 

%% channel rejection?

[~, rejelecs_prob] = pop_rejchan( EEG, 'elec', 1:size(EEG.data,1), ...
    'threshold', 5, 'measure', 'prob', 'norm', 'on');

[~, rejelecs_kurt] = pop_rejchan( EEG, 'elec', 1:size(EEG.data,1), ...
    'threshold', 5, 'measure', 'kurt', 'norm', 'on');

rejelecs_all = union(rejelecs_prob, rejelecs_kurt);

EEG = eeg_checkset( EEG );
EEG = pop_select( EEG, 'nochannel', rejelecs_all);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');

% re-reference head electrodes to average of head electrodes
EEG = eeg_checkset( EEG );
%EEG = pop_reref( EEG, [],'exclude', nonhead_chans);
EEG = pop_reref( EEG, []);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');

%% redo automatic marks with new settings

% clear exisiting rejection marks (if necessary) and re-run all methods
% with new settings
EEG = redo_rejects(EEG, head_noneye_chans_str, head_chans_str, uvthresh, sdthresh);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

eeglab redraw

print_typefreqs(EEG);

%% save current marks

% independently save the current marks for this dataset
[EEG, rejstruct] = save_rejects(EEG);

% save in a structure that keeps each experiment separate
EEG.exprej_struct = save_exprejects(EEG);

[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);