function EEG = import_rawcnttrials_eeglab_inline(data, cnt_srate)

%pull out the set name automatically, get the data, and sampling rate
cnt_data=permute(data,[3 1 2]);

%import data, must be channels x timepts x trials
EEG = pop_importdata('data',cnt_data,'dataformat','array',...
    'nbchan',size(cnt_data,1),...
    'srate',cnt_srate);
EEG = eeg_checkset( EEG );

%add channel locations to remaining channels
EEG=pop_chanedit(EEG, 'lookup','/active_projects/matlab_common/standard_1005.elc',...
    'load',{'/active_projects/matlab_common/61chans_ns.ced' 'filetype' 'autodetect'},...
    'changefield',{47 'type' ''},'changefield',{47 'datachan' 1});
EEG = eeg_checkset( EEG );