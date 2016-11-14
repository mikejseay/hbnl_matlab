function [trial_mat, artifact_vec] = check_artifact(data, trial_mat, ...
		chan_vec, thresh, filter)
    
    if nargin < 5
        filter = false;
    end
    
    if filter
        %import dataR (samps x chans x trials) into eeglab as an EEG structure
        EEG = import_hdftrials_eeglab_inline('12345678901234567', data, 256);
        EEG = eeg_checkset( EEG );

        % filter the data?
        EEG = pop_eegfiltnew(EEG, 20, []);
    end
    
	[n_samps, n_trials, n_chans] = size(data);
	if isempty(chan_vec)
		chan_vec = 1:n_chans;
	end	
	artifact_vec = true(n_trials, 1);
	for m = 1:n_trials
		X = squeeze(data(:, m, chan_vec));
		maxX = max(X);
		minX = min(X);
		diffX = maxX - minX;
		if any(diffX > thresh)
			artifact_vec(m) = false;
		end
	end
	trial_mat = trial_mat & repmat(artifact_vec, [1, size(trial_mat, 2)]);
	artifact_vec = ~artifact_vec;
