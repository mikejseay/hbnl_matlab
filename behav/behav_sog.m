function [avgbet_po,stdbet_po,crit,crit_po,avgbet,resp_mat, ...
 po_mat,avgbet_po2,stdbet_po2] = behav_sog(etable)

% interpet the event sequence in the single-outcome gambling task

% available outs so far
% avgbet_po,stdbet_po,crit,crit_po,avgbet,resp_mat
% po_mat,avgbet_po2,stdbet_po2

n_events=size(etable,1);
num_prev_outcomes_resp=2;
num_prev_outcomes_outcome=2;

outcome_types=[20,40,60,80];
resp_codes=[1,8];
resp_betamts=[10,50];
amount_map=containers.Map(outcome_types,[10,50,-10,-50]);
combamount_map=containers.Map([-100 -60 -40 -20 0 20 40 60 100],[1 1 2 2 3 4 4 5 5]);
resp_map=containers.Map(resp_codes,resp_betamts);
resp_catmap=containers.Map(resp_codes,[1,2]);

resp_mat=zeros(1,2);
prevamount_mat=zeros(1,4);
prevamount2_mat=zeros(1,5);
po_mat=zeros(1,5);
%prevnum_mat=zeros(1,4);
crit_po=zeros(1,4);

prevamount_dum=0;
prevamount2_dum=0;
%for each event
for ev=1:n_events
    
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
                
                %record the amount and count it
                prevamount_dum=prevamount_dum+1;
                prevamount_mat(prevamount_dum, etable.type(ev - checkprev_ev) == outcome_types) = ...
                    resp_map(etable.response_code(ev));
                prevev_tocheck=ev - checkprev_ev;
                             
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
                po_mat(combamount_map(combamount))= ...
                    po_mat(combamount_map(combamount)) + 1;
                             
            end
            
        end
        
    end
    
end

prevamount_mat(prevamount_mat==0)=NaN;
avgbet_po=nanmean(prevamount_mat,1);
stdbet_po=nanstd(prevamount_mat,1);
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
    crit_po(prevoutcome)=-0.5*( norminv(tempprobs(1),0,1) - norminv(tempprobs(2),0,1) );
end

avgbet=(10*resp_mat(1)+50*resp_mat(2))/sum(resp_mat);
prevamount2_mat(prevamount2_mat==0)=NaN;
avgbet_po2=nanmean(prevamount2_mat,1);
stdbet_po2=nanstd(prevamount2_mat,1);

end