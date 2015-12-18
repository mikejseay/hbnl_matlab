function [trial_mat, artifact_vec] = check_artifact(data, trial_mat, ...
		chan_vec, thresh);
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
