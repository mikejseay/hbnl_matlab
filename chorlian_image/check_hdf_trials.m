function [Y, ev_nums] = check_hdf_trials(h1_struct)
	n_trials = h1_struct.experiment_struct.n_trials;
	m_trials_present = numel(h1_struct.trial_struct.artf_present);
    if m_trials_present < n_trials
        Y = [];
        return;
    end
    
	Y = zeros(n_trials, 1);
    for m = 1:n_trials
        if h1_struct.trial_struct.artf_present(m) ~= 1 ...
                && h1_struct.trial_struct.correct(m) ~= 0
            Y(m) = h1_struct.trial_struct.case_num(m);
        end
    end
    ev_nums = find(h1_struct.event_struct.event_id ~= 0)';
end