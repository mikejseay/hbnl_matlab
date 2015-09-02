function [EventLogDS,CNT_Event]=Extract_Events_CNT(file_path)

Resamp_Rate = 256;
Trial_Init_Type=90;

EventLogHeader = {'TrialNum','EventNum','StimType','Keyboard','KeypadAccept', ...
    'AcceptEv1','TimeOffsetOrig','Type','Code','Latency','EpochEvent', ...
    'Accept','Accuracy','TimeOffsetNewRate','TimeLag'};
    
fprintf('Currently processing the file >>  %s [%s]\n',file_path,datestr(now))
S1 = loadcnt(file_path);        % This program is from EEGLAB

CNT_Event = S1.event;
SampRate = S1.header.rate;
SampInt = 1000/SampRate;
Nevent = length(CNT_Event);

Nparameter = length(fieldnames(CNT_Event));	% There are 11 parameters in all CNT files
Event_Mat = NaN(Nevent,Nparameter);

% From the data
for Nev = 1:Nevent
    Event_Mat(Nev,1) = CNT_Event(1,Nev).stimtype;
    Event_Mat(Nev,2) = CNT_Event(1,Nev).keyboard;
    Event_Mat(Nev,3) = CNT_Event(1,Nev).keypad_accept;
    Event_Mat(Nev,4) = CNT_Event(1,Nev).accept_ev1;
    Event_Mat(Nev,5) = CNT_Event(1,Nev).offset;
    Event_Mat(Nev,6) = CNT_Event(1,Nev).type;
    Event_Mat(Nev,7) = CNT_Event(1,Nev).code;
    Event_Mat(Nev,8) = CNT_Event(1,Nev).latency;
    Event_Mat(Nev,9) = CNT_Event(1,Nev).epochevent;
    Event_Mat(Nev,10) = CNT_Event(1,Nev).accept;
    Event_Mat(Nev,11) = CNT_Event(1,Nev).accuracy;
end

% Prefix a column "Event_Num"
Event_Num = (1:Nevent)';
Event_Mat = horzcat(Event_Num, Event_Mat);

% Prefix a column  as the first column having "Trial_Num" based on CHOICE stimuli
Stim_index = find(Event_Mat(:,2)==Trial_Init_Type);
Trial_Num = NaN(size(Event_Mat,1),1);		% Create an NaN matrix
for ST1 = 1:length(Stim_index)
    Trial_Num(Stim_index(ST1)) = ST1;
end

for EI1 = 2:size(Event_Mat,1)		% Fill the remaining NaN with previous number
    if isnan(Trial_Num(EI1))==1
        Trial_Num(EI1) = Trial_Num(EI1-1);
    end
end
Event_Mat = horzcat(Trial_Num,Event_Mat);	% Prefix Trial_Num as the first column

% Adding a TF_Offset_New
TF_Offset_NewRate = round(Event_Mat(:,7).*(Resamp_Rate/SampRate));
Event_Mat = horzcat(Event_Mat,TF_Offset_NewRate);

% Computing and adding a column the Timelag vector at the end
TF_Lag = zeros(Nevent,1);
for N3 = 2:length(TF_Lag)
    TF_Lag(N3) = (Event_Mat(N3,7))-(Event_Mat((N3-1),7));
end
Time_Lag = TF_Lag .* SampInt;
Event_Mat = horzcat(Event_Mat,Time_Lag);

% Remove 'NaN' rows and correct the second column (due to recording bugs related to "late" start,
% e.g., trials starting with 'button-press' or 'outcome stimulus' as in
% ern_8_b1_20121011_32.cnt)
Event_Mat = Event_Mat(~isnan(Event_Mat(:,1)),:);
Event_Mat(:,2) = 1:size(Event_Mat,1);

% Recomputing 'accuracy' column
Event_Mat(:,13)=0;      % First, make all values to '0'

Outcome_index = find((Event_Mat(:,3)==80) | (Event_Mat(:,3)==60) | (Event_Mat(:,3)==40) | (Event_Mat(:,3)==20));

Choice_index = zeros(length(Outcome_index),1);      % Choice index for trials with feedback
for Noi = 1:length(Outcome_index)
    Choice_index(Noi) = find((Event_Mat(:,1)==(Event_Mat(Outcome_index(Noi),1)) & (Event_Mat(:,3)==90)));
end

for Nidx = 1:length(Outcome_index)        % This Choice_index fixed the problem of multiple button presses within a trial
    if ((Event_Mat((Choice_index(Nidx)+1),3)==0) && ((Event_Mat((Choice_index(Nidx)+1),5)==1) || (Event_Mat((Choice_index(Nidx)+1),5)==8)) && (Event_Mat((Choice_index(Nidx)+1),15)>=100) && (Event_Mat((Choice_index(Nidx)+1),15)<=1000))
        Event_Mat(Outcome_index(Nidx),13) = 1;
    else
        Event_Mat(Outcome_index(Nidx),13) = -1;
    end
end
EventLogDS = mat2dataset(Event_Mat,'VarNames',EventLogHeader);