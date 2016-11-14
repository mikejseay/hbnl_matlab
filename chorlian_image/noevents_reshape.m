function [trial_index, data] = noevents_reshape(Xdat, n_samps) %, inrate, outrate)

[total_samps, n_chans] = size(Xdat);
n_trials = floor(total_samps ./ n_samps);

trial_index = zeros(n_trials, 1);
for m = 1:n_trials
    trial_index(m) = 1 + n_samps.*(m-1);
end

data = zeros(n_samps, n_trials, n_chans);
for m = 1:n_trials
    data(:, m, :) = ...
        Xdat(trial_index(m):trial_index(m) + n_samps - 1, :);
end