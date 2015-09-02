function [trial_mat,typedescriptors,respRT] = behavinds_sog(etable)

% interpet the event sequence in the single-outcome gambling task

% available outs so far
% avgbet_prevoutcome,stdbet_prevoutcome,crit,crit_prevoutcome,avgbet,resp_mat
% prevoutcome_mat, avgbet2_prevoutcome,stdbet2_prevoutcome

% response: 1=bet of 10, 8=bet of 50

n_events=size(etable,1);
num_prev_outcomes_resp=1;
num_prev_outcomes_outcome=1;

outcome_types=[20,40,60,80];
resp_codes=[1,8];
resp_betamts=[10,50];
amount_map=containers.Map(outcome_types,[10,50,-10,-50]);
combamount_map=containers.Map([-100 -60 -40 -20 0 20 40 60 100],[1 1 2 2 3 4 4 5 5]);
resp_map=containers.Map(resp_codes,resp_betamts);
resp_catmap=containers.Map(resp_codes,[1,2]);
outcome_catmap=containers.Map(outcome_types,[1 1 2 2]);

resp_mat=zeros(1,2);
prevamount_mat=zeros(1,4);
prevamount2_mat=zeros(1,5);
prevoutcome_mat=zeros(1,5);
%prevnum_mat=zeros(1,4);
crit_prevoutcome=zeros(1,4);
respRT_mat=zeros(1,4);
respRT=zeros(3,4);

n_trials=sum(etable.type~=0);
n_types=9;
typedescriptors={'Bet 10 Post-Loss','Bet 10 Post-Win','Bet 50 Post-Loss','Bet 50 Post-Win',...
    '-80 Trend','-30 Trend','0 Trend','30 Trend','80 Trend'};
type_indcombs={[1 1],[1 2],[2 1],[2 2]};
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
    
    %if a valid response
    if ismember(etable.response_code(ev),resp_codes) && etable.errant_response(ev)==0 && ...
            ( etable.rt(ev) > 99 && etable.rt(ev) < 1001)
        
        
        %count the response choice
        resp_mat(resp_catmap(etable.response_code(ev))) = ...
            resp_mat(resp_catmap(etable.response_code(ev))) + 1;
        
        %for some number of previous *outcomes*
        prevoutcome_dum=0;
        checkprev_ev=0;
        prevev_tocheck=0;
        while prevoutcome_dum<num_prev_outcomes_resp
            checkprev_ev=checkprev_ev+1;
            
            if prevoutcome_dum==0
            %if event index < 1, break the while loop
            if ev - checkprev_ev < 1
                break
            end
            
            %check to see if the previous event was an outcome
            if ismember(etable.type(ev - checkprev_ev),outcome_types)
                prevoutcome_dum=prevoutcome_dum+1;
                respRT_dum=respRT_dum+1;
                
                %record the amount and count it
                prevamount_dum=prevamount_dum+1;
                prevamount_mat(prevamount_dum, etable.type(ev - checkprev_ev) == outcome_types) = ...
                    resp_map(etable.response_code(ev));
                prevev_tocheck=ev - checkprev_ev;
                
                %index!
                type_inds=[resp_catmap(etable.response_code(ev)) ...
                    outcome_catmap(etable.type(ev - checkprev_ev))];
                resp_ind=find(cellfun(@all,cellfun(@(x) x==type_inds,...
                    type_indcombs,'UniformOutput',false)));
                trial_mat(trial_dum,resp_ind) = true;
                
                %record the RT
                respRT_mat(respRT_dum,resp_ind)=etable.rt(ev);
                             
            end
            elseif prevoutcome_dum==1
            %if event index < 1, break the while loop
            if ev - checkprev_ev < 1
                break
            end
            
            %check to see if the previous event was an outcome
            if ismember(etable.type(ev - checkprev_ev),outcome_types)
                prevoutcome_dum=prevoutcome_dum+1;
                
                %calculate the combined previous outcomes amount, and
                %record the bet following it
                prevamount2_dum=prevamount2_dum+1;
                combamount = amount_map(etable.type(prevev_tocheck)) + ...
                    amount_map(etable.type(ev - checkprev_ev));
                prevamount2_mat(prevamount2_dum, combamount_map(combamount)) = ...
                    resp_map(etable.response_code(ev));
            
            end
            end
            
        end
        
    end
    
    %if an outcome
    if ismember(etable.type(ev),outcome_types)
        
        %for some number of previous *outcomes*
        prevoutcome_dum=0;
        checkprev_ev=0;
        while prevoutcome_dum<num_prev_outcomes_outcome
            checkprev_ev=checkprev_ev+1;
            
            %if event index < 1, break the while loop
            if ev - checkprev_ev < 1
                break
            end
            
            %check to see if the previous event was an outcome
            if ismember(etable.type(ev - checkprev_ev),outcome_types)
                prevoutcome_dum=prevoutcome_dum+1;
                
                %record the combined amount and categorize it
                combamount = amount_map(etable.type(ev)) + ...
                    amount_map(etable.type(ev - checkprev_ev));
                prevoutcome_mat(combamount_map(combamount))= ...
                    prevoutcome_mat(combamount_map(combamount)) + 1;
                
                %index!
                trial_mat(trial_dum, ...
                    combamount_map(combamount)+4) = ...
                    true;
                             
            end
            
        end
        
    end
    
end

prevamount_mat(prevamount_mat==0)=NaN;
avgbet_prevoutcome=nanmean(prevamount_mat,1);
stdbet_prevoutcome=nanstd(prevamount_mat,1);
respprob=resp_mat/(sum(resp_mat));
if any(respprob==0 | respprob==1)
    respprob(respprob==0)=1/(2*sum(resp_mat));
    respprob(respprob==1)=(2*sum(resp_mat)-1)/(2*sum(resp_mat));
end
crit=-0.5*( norminv(respprob(1),0,1) - norminv(respprob(2),0,1) );

for prevoutcome=1:size(prevamount_mat,2)
    tempprobs=zeros(1,2);
    for resp=1:2
        tempprobs(resp)=sum(prevamount_mat(:,prevoutcome)==resp_betamts(resp)) / ...
            sum(~isnan(prevamount_mat(:,prevoutcome)));
    end
    if any(tempprobs==0 | tempprobs==1)
        tempprobs(tempprobs==0)=1/(2*sum(~isnan(prevamount_mat(:,prevoutcome))));
        tempprobs(tempprobs==1)=(2*sum(~isnan(prevamount_mat(:,prevoutcome)))-1) / ...
            (2*sum(~isnan(prevamount_mat(:,prevoutcome))));
    end
    crit_prevoutcome(prevoutcome)=-0.5*( norminv(tempprobs(1),0,1) - norminv(tempprobs(2),0,1) );
end

avgbet=(10*resp_mat(1)+50*resp_mat(2))/sum(resp_mat);
prevamount2_mat(prevamount2_mat==0)=NaN;
avgbet2_prevoutcome=nanmean(prevamount2_mat,1);
stdbet2_prevoutcome=nanstd(prevamount2_mat,1);

respRT_mat(respRT_mat==0)=NaN;
respRT(1,:)=nanmean(respRT_mat,1);
respRT(2,:)=nanmedian(respRT_mat,1);
respRT(3,:)=nanvar(respRT_mat,1);

end