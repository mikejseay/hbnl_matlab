function EEG = redo_rejects(EEG, threshchans, probkurtchans, uvthresh, sdthresh)

threshchan_inds = convert_chanlabels(EEG.chanlocs, threshchans);
probkurtchan_inds = convert_chanlabels(EEG.chanlocs, probkurtchans);

% clear existing reject marks
EEG = clear_rejects(EEG);

EEG = pop_eegthresh(EEG, 1, threshchan_inds, -uvthresh, uvthresh, EEG.xmin, EEG.xmax, 1, 0);
EEG = pop_jointprob(EEG, 1, probkurtchan_inds, sdthresh, sdthresh, 1, 0);
EEG = pop_rejkurt(EEG, 1, probkurtchan_inds, sdthresh, sdthresh, 1, 0);

end