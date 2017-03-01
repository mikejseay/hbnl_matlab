function trial_mat_out = check_artifact_eeglab(data_file, data, h1_struct, ...
    trial_mat_in, sdthresh, uvthresh, chans)

%import dataR (samps x chans x trials) into eeglab as an EEG structure
EEG = import_hdftrials_eeglab_inline(data_file, data, h1_struct);
EEG = eeg_checkset( EEG );

% filter the data?
%EEG = pop_eegfiltnew(EEG, 30, []);

%run the automatic epoch rejection
[EEG, rmepochs] = pop_autorej(EEG, 'threshold', uvthresh, 'electrodes', chans, ...
    'startprob', sdthresh, 'maxrej', 20, 'nogui', 'on');

%use the output to mask the trial indexing scheme
rmepochs=sort(rmepochs);
trial_mat_out=trial_mat_in;
trial_mat_out(rmepochs,:)=false;

end