function [trial_mat, artifact_vec] = check_artifact_baselined(data, trial_mat, ...
		chan_vec, thresh, baseline_pts)
	[n_samps, n_trials, n_chans] = size(data);
	if isempty(chan_vec)
		chan_vec = 1:n_chans;
	end	
	artifact_vec = true(n_trials, 1);
    % should apply lopass filter here
	for m = 1:n_trials
		X = squeeze(data(:, m, chan_vec));
        XX = bsxfun(@minus, X, mean(X(baseline_pts, :)));
		absmaxXX = max(abs(XX));
		if any(absmaxXX  > thresh)
			artifact_vec(m) = false;
		end
	end
	trial_mat = trial_mat & repmat(artifact_vec, [1, size(trial_mat, 2)]);
	artifact_vec = ~artifact_vec;
