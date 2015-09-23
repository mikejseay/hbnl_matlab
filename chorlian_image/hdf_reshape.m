function [trial_index, n_samps, data] = hdf_reshape(Xdat, h1_struct, pstruct)
	[total_samps, n_chans] = size(Xdat);
	n_trials = h1_struct.experiment_struct.n_trials;
	samp_rate = h1_struct.experiment_struct.rate;
	trial_index = zeros(n_trials, 1);
	if pstruct.rate ~= samp_rate % resample if necessary 
		Xdat = resample(Xdat, pstruct.rate, samp_rate);
	end
	for m = 1:n_trials
		trial_index(m) = ...
			round(h1_struct.trial_struct.time_offset(m) * pstruct.rate);
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
	if trial_index(n_trials) + n_samps - 1 > total_samps
		trial_index = trial_index(1:(end - 1));
		n_trials = n_trials - 1;
	end	
	if nargout < 3
		return
    end
    % factor in baseline length
    baseline_pts = round(pstruct.prestim_ms / 1000 * pstruct.rate);
    % if very long, clip it back a bit
    if n_samps > pstruct.n_samps + 10
        n_samps = pstruct.n_samps + 10;
    end
    % if want to lengthen baseline onto end
    if pstruct.lengthen
        n_samps = n_samps + baseline_pts;
    end    
    %
	data = zeros(n_samps, n_trials, n_chans);
    %check for complete data from first trial
    if (trial_index(1) - baseline_pts) < 1
        start_m=2;
    else
        start_m=1;
    end
	for m = start_m:n_trials
		data(:, m, :) = ...
			Xdat((trial_index(m) - baseline_pts):(trial_index(m) + n_samps - baseline_pts - 1), :);
    end

	
