function [trial_type, ev_nums] = check_hdf_trials2(h1_struct)
	event_types = h1_struct.event_struct.event_id;
    n_events = length(event_types);
    trials = event_types(event_types~=0);
    n_trials = length(trials);
    
	m_trials_present = numel(h1_struct.trial_struct.artf_present);
    if m_trials_present < n_trials
        trial_type = [];
        return;
    end
    
	trial_type = zeros(n_trials, 1);
    ev_nums = zeros(n_events, 1);
	for m = 1:n_trials
		if h1_struct.trial_struct.artf_present(m) ~= 1 ...
				&& h1_struct.trial_struct.correct(m) ~= 0
			trial_type(m) = h1_struct.trial_struct.case_num(m);
		end
	end
end