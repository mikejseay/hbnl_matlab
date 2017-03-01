function ev_log = make_evlog(n_events, prop_true, n_consec)

n_wins = floor(n_events * prop_true);
%n_losses = n_events - n_wins;

win_inds = sort(randperm(n_events, n_wins));

ev_log = false(n_events, 1);
ev_log(win_inds) = true;

ev_log = fix_consecutive_boolean(ev_log, n_consec);

end