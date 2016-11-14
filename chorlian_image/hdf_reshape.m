function [trial_index, n_samps, data, trial_type_out, ev_num_out] = ...
    hdf_reshape(Xdat, h1_struct, pstruct, trial_type_in, ev_num_in)
	[total_samps, n_chans] = size(Xdat);
	n_trials = h1_struct.experiment_struct.n_trials;
    trial_type_out = trial_type_in;
    ev_num_out = ev_num_in;
    rate = h1_struct.experiment_struct.rate;
	
    trial_index = round(h1_struct.trial_struct.time_offset * rate);
    
	n_samps = round( sum(abs( pstruct.epoch_lims )) ./ 1000 * rate );
    
    % factor in baseline length
    baseline_pts = round(-pstruct.epoch_lims(1) ./ 1000 * rate);
        
	% check for complete data from last trials
    while trial_index(n_trials) + n_samps - baseline_pts - 1 > total_samps
		trial_index = trial_index(1:(end - 1));
        trial_type_out = trial_type_out(1:(end - 1));
        ev_num_out = ev_num_out(1:(end - 1));
		n_trials = n_trials - 1;
    end
    % check for complete data from first trials
    while trial_index(1) - baseline_pts < 1
		trial_index = trial_index(2:end);
        trial_type_out = trial_type_out(2:end);
        ev_num_out = ev_num_out(2:end);
		n_trials = n_trials - 1;
    end
    
    % preallocate the matrix
	data = zeros(n_samps, n_trials, n_chans);
    
    for m = 1:n_trials
		data(:, m, :) = ...
			Xdat((trial_index(m) - baseline_pts):(trial_index(m) + n_samps - baseline_pts - 1), :);
    end
end