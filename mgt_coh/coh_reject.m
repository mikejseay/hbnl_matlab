%% data rejection step

s_logic=false(imp.s_valid,3);

%bit to pick out bad data based on high coherence values
% catches channel pairs that are bridged, using a
% simple threshold

%%
if false
%if true
high_coh_vals=zeros(imp.s_valid,imp.maxpairs);
for chosen_s=1:imp.s_valid
    for pair=1:imp.maxpairs
        chosen_s_data=cohdata(:,6:15,:,pair,chosen_s);
        high_coh_vals(chosen_s,pair)=numel(chosen_s_data(chosen_s_data>0.9999)); %/numel(chosen_s_data);
    end
    if sum(high_coh_vals(chosen_s,:),2)>50*imp.maxpairs
        s_logic(chosen_s,1)=false;
    else
        s_logic(chosen_s,1)=true;
    end
end        
high_coh_vals_by_s=sum(high_coh_vals,2);


% bit to pick out bad data based on highly correlated coherence values
% -catches channel pairs which share a common bridged electrode (i.e.
% although FZ-O1 and FZ-OZ both have plausible coherences (<<.99), they
% are near-identical (near-perfectly correlated) because O1 and OZ are 
% bridged to each other, using a simple threshold

s_highrhocount=zeros(imp.s_valid,1);
for chosen_s=1:imp.s_valid
    chosen_s_data=squeeze(cohdata(:,10,1,:,chosen_s));    
    %
    rho=tril(corr(chosen_s_data),-1);
    %
    s_highrhocount(chosen_s)=numel(rho(rho>.995));
    if  s_highrhocount(chosen_s)>round(imp.maxpairs/2);
        s_logic(chosen_s,2)=false;
    else
        s_logic(chosen_s,2)=true;
    end
end
end
%%

%bit to pick out extreme/implausible ITC values??

%s_highitccount=zeros(imp.s_valid,1);
%for chosen_s=1:imp.s_valid
%    chosen_s_data=squeeze(mean(itcdata(:,:,6:20,1,chosen_s),1));
%    s_highitccount(chosen_s) = numel(chosen_s_data(chosen_s_data > .9));
%    s_itcmat(chosen_s)=chosen_s_data;
%end

%bit to pick out extreme/implausible ERO values??

%bit to pick out extreme/implausible ERP values??

%bit to exclude subjects which have insufficient data (i.e. not enough
%trials per condition)
%trial_thresh=ceil(mean(n_trials_all,2)-std(n_trials_all,0,2));
if true
trials_necessary=15;
trial_thresh=trials_necessary*ones(imp.maxconds,1);
for chosen_s=1:imp.s_valid
    s_logic(chosen_s,3)=all(n_trials_all(:,chosen_s) > trial_thresh);
end
end

% combine rejection logic
%s_inds=false(imp.s_valid,1);
%for chosen_s=1:imp.s_valid
%    if any(s_logic(chosen_s,:))
%        s_inds(chosen_s)=1;
%    end
%end
s_inds=any(s_logic,2);

%custom rejection
if false
    custom_bad_s=[23    47    48    49    66    77    86    88    92    95    99 ...
        106   110   116   121   124   134   137   138   143   145   146   153   158 ...
        162   165   167   170   173   174   176   180   181   186   189   198   207 ...
        216   224   225   227];
    custom_dummy=0;
    for chosen_s=1:imp.s_valid
        if s_inds(chosen_s)
            custom_dummy=custom_dummy+1;
            if ismember(custom_dummy,custom_bad_s)
                s_inds(chosen_s)=false;
            end
        end
    end
end

fprintf('%d percent of subjects remaining (%d / %d)\n',round(sum(s_inds)/imp.s_valid*100),sum(s_inds),imp.s_valid)

clear high_coh_vals chosen_s_data rho trial_thresh chosen_s pair custom_bad_s custom_dummy s_logic trials_necessary