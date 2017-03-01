function trial_mat_out = check_artifact_eeglab_noevents(data, cnt_srate, ...
    trial_mat_in, sdthresh, uvthresh, chans)

%import dataR (samps x chans x trials) into eeglab as an EEG structure
EEG = import_rawcnttrials_eeglab_inline(data, cnt_srate);
EEG = eeg_checkset( EEG );

%run the automatic epoch rejection
[EEG, rmepochs] = pop_autorej(EEG, 'threshold', uvthresh, 'electrodes', chans, ...
    'startprob', sdthresh, 'maxrej', 15, 'nogui', 'on');

%use the output to mask the trial indexing scheme
rmepochs=sort(rmepochs);
trial_mat_out=trial_mat_in;
trial_mat_out(rmepochs,:)=false;

end