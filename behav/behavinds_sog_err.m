function [trial_mat,typedescriptors,respRT] = behavinds_sog_err(etable)

% interpet the event sequence in the single-outcome gambling task

% response: 1=bet of 10, 8=bet of 50

n_events=size(etable,1);
num_prev_outcomes_resp=1;

outcome_types=[20,40,60,80];
resp_codes=[1,8];
resp_betamts=[10,50];
amount_map=containers.Map(outcome_types,[10,50,-10,-50]);
resp_map=containers.Map(resp_codes,resp_betamts);
resp_catmap=containers.Map(resp_codes,[1,2]);
outcome_catmap=containers.Map(outcome_types,[1 1 2 2]);

resp_mat=zeros(1,2);
respRT_mat=zeros(1,4);
respRT=zeros(3,4);

n_trials=sum(etable.type~=0);
n_types=2;
typedescriptors={'Gain','Loss'};
trial_mat=false(n_trials,n_types);

prevamount_dum=0;
prevamount2_dum=0;
trial_dum=0;
respRT_dum=0;
%for each event
for ev=1:n_events
    
    if etable.type(ev)~=0
        trial_dum=trial_dum+1;
    end
    
    %if a valid response (RT = [100,1000])
    if ismember(etable.response_code(ev),resp_codes) && etable.errant_response(ev)==0 && ...
            ( etable.rt(ev) >= 100 && etable.rt(ev) <= 1000)
        
        
        %count the response choice
        resp_mat(resp_catmap(etable.response_code(ev))) = ...
            resp_mat(resp_catmap(etable.response_code(ev))) + 1;
        
    end
    
    %if an outcome
    if ismember(etable.type(ev),outcome_types)
        
        %if the previous response was valid
        
        
        %index it
        trial_mat(trial_dum, outcome_catmap( etable.type(ev) )) = true;
    end
    
end

respRT_mat(respRT_mat==0)=NaN;
respRT(1,:)=nanmean(respRT_mat,1);
respRT(2,:)=nanmedian(respRT_mat,1);
respRT(3,:)=nanvar(respRT_mat,1);

end