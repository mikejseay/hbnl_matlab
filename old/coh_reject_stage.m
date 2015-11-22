function s_inds = coh_reject_stage (imp, n_trials_all, trials_necessary)
% subject rejection step

if nargin < 3
    trials_necessary=15;
end

s_logic=false(imp.s_valid,3);

trial_thresh=trials_necessary*ones(imp.maxconds,1);
for chosen_s=1:imp.s_valid
    s_logic(chosen_s,3)=all(n_trials_all(:,chosen_s) > trial_thresh);
end

s_inds=any(s_logic,2);

fprintf('%d percent of subjects remaining (%d / %d)\n', ...
    round(sum(s_inds)/imp.s_valid*100),sum(s_inds),imp.s_valid)

end