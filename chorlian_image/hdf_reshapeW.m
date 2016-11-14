function [trial_index, n_samps, data] = hdf_reshapeW(Wdat, h1_struct, pstruct)
	[total_samps, n_chans, n_scales] = size(Wdat);
	n_trials = h1_struct.experiment_struct.n_trials;
    rate = h1_struct.experiment_struct.rate;
    
    trial_index = round(h1_struct.trial_struct.time_offset * rate);
    
    n_samps = round( sum(abs( pstruct.epoch_lims )) ./ 1000 * rate );
    
    % factor in baseline length
    baseline_pts = round(-pstruct.epoch_lims(1) ./ 1000 * rate);
        
	% check for complete data from ending trials
    while trial_index(n_trials) + n_samps - baseline_pts - 1 > total_samps
		trial_index = trial_index(1:(end - 1));
		n_trials = n_trials - 1;
    end
    % check for complete data from beginning trials
    while trial_index(1) - baseline_pts < 1
		trial_index = trial_index(2:end);
		n_trials = n_trials - 1;
    end
    
    %preallocate
    data = zeros(n_samps, n_trials, n_chans, n_scales);
    
    for m = 1:n_trials
		data(:, m, :, :) = ...
			Wdat( (trial_index(m) - baseline_pts):(trial_index(m) + n_samps - baseline_pts - 1), :, :);
    end
    
end