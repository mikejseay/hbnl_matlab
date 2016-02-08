function s_inds = coh_reject (imp, n_trials_all, trials_necessary, custom_rej)
% subject rejection step

if nargin < 4
    custom_rej=[];
end
if nargin < 3
    trials_necessary=15;
end

s_logic=false(imp.s_valid,2);

trial_thresh=trials_necessary*ones(imp.maxconds,1);
for chosen_s=1:imp.s_valid
    
    %trials
    s_logic(chosen_s,1)=all(n_trials_all(:,chosen_s) > trial_thresh);
    
    %custom
    if ismember(chosen_s,custom_rej)
        s_logic(chosen_s,2)=false;
    else
        s_logic(chosen_s,2)=true;
    end    
end

s_inds=all(s_logic,2);

end