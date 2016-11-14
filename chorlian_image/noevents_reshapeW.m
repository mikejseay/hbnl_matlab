function [trial_index, data] = noevents_reshapeW(Wdat, n_samps)  
[total_samps, n_chans, n_scales] = size(Wdat);
n_trials = floor(total_samps ./ n_samps);

trial_index = zeros(n_trials, 1);
for m = 1:n_trials
    trial_index(m) = 1 + n_samps.*(m-1);
end

data = zeros(n_samps, n_trials, n_chans, n_scales);
for m = 1:n_trials
    data(:, m, :, :) = ...
        Wdat(trial_index(m):trial_index(m) + n_samps - 1, :, :);
end