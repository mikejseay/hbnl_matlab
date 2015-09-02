function [trial_index, n_samps, data] = hdf_reshapeW(Wdat, h1_struct, rate, prestim_ms)
	[total_samps, n_chans, n_scales] = size(Wdat);
	n_trials = h1_struct.experiment_struct.n_trials;
	samp_rate = h1_struct.experiment_struct.rate;
	trial_index = zeros(n_trials, 1);
    
	for m = 1:n_trials
		trial_index(m) = ...
			 round(h1_struct.trial_struct.time_offset(m) * rate);
    end
	[n_samps, idx] = min(diff(trial_index));
	while n_samps == 0
		n_trials = idx - 1;
		[n_samps, idx] = min(diff(trial_index(1:n_trials)));
	end
	if n_trials == 0
		data = [];
		trial_index = [];
		return	
    end
	% check for complete data from last trial
	if trial_index(n_trials) + n_samps - 1 > total_samps
		trial_index = trial_index(1:(end - 1));
		n_trials = n_trials - 1;
	end
	if nargout < 3
		return
    end
    % factor in baseline length
    baseline_pts = round(prestim_ms / 1000 * rate);
    %n_samps = n_samps + baseline_pts;
    %
	data = zeros(n_samps, n_trials, n_chans, n_scales);
    % check for complete data from first trial
    if (trial_index(1) - baseline_pts) < 1
        start_m=2;
    else
        start_m=1;
    end
	for m = start_m:n_trials
		data(:, m, :, :) = ...
			Wdat( (trial_index(m) - baseline_pts) :(trial_index(m) + n_samps - baseline_pts - 1), :, :);
    end