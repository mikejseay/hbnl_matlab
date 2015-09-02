function [TRBvarsDS]=Compute_TRBvars_Kam(file_list,dir_path,extension)

%% This program computes TRB variables directly from the CNT files (by creating an 'event-log' from each) [@created by Chella Kamarajan in Sep'2012]

%% Clock Starts
tic
%% PART I : Signal Processing Options

if nargin>1

% Get the list of CNT-files
GNG_CNT_List = getAllFiles2(dir_path, extension,1);

% Select the list of 'Phase4' subjects (using 'regexp')
nCNTall = numel(GNG_CNT_List);
GNG_CNT_List_Ph4 = cell(size(GNG_CNT_List));
for iCNT = 1:nCNTall
    GNG_CNT_List_Ph4{iCNT} =  regexp(GNG_CNT_List{iCNT}, 'ern_\d{1}_\w{1}\d{1}_\d{8}_32.cnt', 'match');
end
GNG_CNT_Ph4_Path_idx = find(~cellfun(@isempty,GNG_CNT_List_Ph4));
GNG_CNT_Ph4_Path = GNG_CNT_List(GNG_CNT_Ph4_Path_idx);

% Save the file-list as a CSV file
export((cell2dataset(GNG_CNT_Ph4_Path,'ReadVarNames',0)),'File','MGT_CNT_Ph4_List_n7603.csv','WriteVarNames',0,'Delimiter',',');

end

GNG_CNT_Ph4_Path=file_list;

%% Create lists of File_ID, Folder_ID, and Subject_ID
Conditions = {'N50','N10','P50','P10'};
Stim_Codes = [80,60,40,20];
Resamp_Rate = 256;

nSubjects = length(GNG_CNT_Ph4_Path);
nConditions = length(Conditions);
FileIDs = cell(nSubjects,1);
SubjectIDs = cell(nSubjects,1);
SessRun = cell(nSubjects,1);
Session = cell(nSubjects,1);

for iSubj = 1:nSubjects
    [~,FileIDs{iSubj},~] = fileparts(GNG_CNT_Ph4_Path{iSubj});
    FileIDs{iSubj} = FileIDs{iSubj}(1:17);
    SubjectIDs{iSubj} = FileIDs{iSubj}(10:17);
    SessRun{iSubj} = FileIDs{iSubj}(7:8);
    Session{iSubj} = FileIDs{iSubj}(7);
end
Task = upper(FileIDs{1}(1:3));

%% TRB Variables
TRBvarsList = {'N_Tot_Resp','N50_nTrl','N10_nTrl','P50_nTrl','P10_nTrl','N50_RT_M','N10_RT_M','P50_RT_M','P10_RT_M','Loss_RT_M','Gain_RT_M','Amt50_RT_M','Amt10_RT_M',...
    'SF_50_LS_trials_1tr','SF_50_LS_trials_2tr','SF_50_LS_trials_3tr','SF_50_LS_trials_4tr','SF_50_LS_trials_5tr','SF_50_LS_trend_2tr','SF_50_LS_trend_3tr','SF_50_LS_trend_4tr','SF_50_LS_trend_5tr'};
nTRBvars = numel(TRBvarsList);
TRBvarsDS = [cell2dataset(cell(nSubjects,3),'VarNames',{'SubjID','SessRun','Session'}),mat2dataset(nan(nSubjects,nTRBvars),'VarNames',TRBvarsList)];

%% Compute EventLog and Extract TRB data
EventLogHeader = {'TrialNum','EventNum','StimType','Keyboard','KeypadAccept','AcceptEv1',...
    'TimeOffsetOrig','Type','Code','Latency','EpochEvent','Accept','Accuracy','TimeOffsetNewRate','TimeLag'};
for iSubj = 1:nSubjects
    
    %% Read the cnt_file from the cell
    cnt_file_each = GNG_CNT_Ph4_Path{iSubj};
    
    %% Load the cnt-file (-> output is a structure)
    fprintf('Currently processing the file >>  %d of %d [%s]\n',iSubj,nSubjects,datestr(now))   % Info at the screen
    S1 = loadcnt(cnt_file_each);        % This program is from EEGLAB
    
    %% Parameters from the Data Structure
    SampRate = S1.header.rate;
    SampInt = 1000/SampRate;
    Nelec = S1.header.nchannels;	% This includes the BLANK channel
    Nevent = length(S1.event);
    
    %% Electrodes Information
    % Electrodes in Original Raw Data
    Elec_Labels_All64 = cell(1,Nelec);
    for Nel = 1:Nelec
        Elec_Labels_All64{1,Nel} = S1.electloc(1,Nel).lab;
    end
    
    % Electrodes for signal processing
    Nelec_NoBlank = [(1:62),64];     % 63rd electrode is blank, which does not go into analyses
    Elec_Labels_NoBlank = Elec_Labels_All64(Nelec_NoBlank);
    Elec_Num_Scalp = [(1:31),(33:62)];      % X=32; BLANK=63, Y=64
    Elec_Labels_Scalp = Elec_Labels_All64(Elec_Num_Scalp);
    
    %% Creating the Event Log
    Nparameter = length(fieldnames(S1.event));	% There are 11 parameters in all CNT files
    Event_Mat = NaN(Nevent,Nparameter);
    
    % From the data
    for Nev = 1:Nevent
        Event_Mat(Nev,1) = S1.event(1,Nev).stimtype;		% This could have been done with a cell structure as well
        Event_Mat(Nev,2) = S1.event(1,Nev).keyboard;
        Event_Mat(Nev,3) = S1.event(1,Nev).keypad_accept;
        Event_Mat(Nev,4) = S1.event(1,Nev).accept_ev1;
        Event_Mat(Nev,5) = S1.event(1,Nev).offset;
        Event_Mat(Nev,6) = S1.event(1,Nev).type;
        Event_Mat(Nev,7) = S1.event(1,Nev).code;
        Event_Mat(Nev,8) = S1.event(1,Nev).latency;
        Event_Mat(Nev,9) = S1.event(1,Nev).epochevent;
        Event_Mat(Nev,10) = S1.event(1,Nev).accept;
        Event_Mat(Nev,11) = S1.event(1,Nev).accuracy;
    end
    
    % Prefix a column "Event_Num"
    Event_Num = (1:Nevent)';
    Event_Mat = horzcat(Event_Num, Event_Mat);
    
    % Prefix a column  as the first column having "Trial_Num" based on CHOICE stimuli
    Stim_index = find(Event_Mat(:,2)==90);
    Trial_Num = NaN(size(Event_Mat,1),1);		% Create an NaN matrix
    for ST1 = 1:length(Stim_index)
        Trial_Num(Stim_index(ST1)) = ST1;
    end
    
    for EI1 = 2:size(Event_Mat,1)		% Fill the remaining NaN with previous number
        if isnan(Trial_Num(EI1))==1
            Trial_Num(EI1) = Trial_Num(EI1-1);
        end
    end
    Event_Mat = horzcat(Trial_Num,Event_Mat);	%#ok<*AGROW> % Prefix the 'Trial_Num' as the first column
    
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
    
    % Remove 'NaN' rows and correct the second column (due to the recording bugs related to "late" start,
    % e.g., trials starting with 'button-press' or 'outcome stimulus' as in
    % ern_8_b1_20121011_32.cnt
    Event_Mat = Event_Mat(~isnan(Event_Mat(:,1)),:);
    Event_Mat(:,2) = 1:size(Event_Mat,1);
    
    % Recomputing 'accuracy' column
    Event_Mat(:,13)=0;      % First, make all values to '0'
    
    
    Outcome_index = find((Event_Mat(:,3)==80) | (Event_Mat(:,3)==60) | (Event_Mat(:,3)==40) | (Event_Mat(:,3)==20));
    %     for Nidx = 1:length(Outcome_index)        % This produced erroneous 'accuracy' computations. So, Choice_index is to be tried
    %         if ((Event_Mat((Outcome_index(Nidx)-1),3)==0) && ((Event_Mat((Outcome_index(Nidx)-1),5)==1) || (Event_Mat((Outcome_index(Nidx)-1),5)==8)) && (Event_Mat((Outcome_index(Nidx)-1),15)>=100) && (Event_Mat((Outcome_index(Nidx)-1),15)<=1000))
    %             Event_Mat(Outcome_index(Nidx),13) = 1;
    %         else
    %             Event_Mat(Outcome_index(Nidx),13) = -1;
    %         end
    %     end
    
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
    %if iSubj==1
        %export(EventLogDS,'file',sprintf('%s_events.csv', FileIDs{iSubj}),'Delimiter',',');
    %end
        
    %% Selection Frequency for Amount-50 when previous trials had a 'losing trend'
    % Create 'OutcomeDS' with only the the outcome events
    OutcomeDS = EventLogDS((EventLogDS.StimType==80 | EventLogDS.StimType==60 | EventLogDS.StimType==40 | EventLogDS.StimType==20),:);
    OutcomeDS.SrlNum = transpose(1:size(OutcomeDS,1));
    OutcomeDS.StimVal = OutcomeDS.StimType;
    OutcomeDS = OutcomeDS(:,[16,1:3,17,4:15]);
    OutcomeDS.StimVal(OutcomeDS.StimType==80) = -50;
    OutcomeDS.StimVal(OutcomeDS.StimType==60) = -10;
    OutcomeDS.StimVal(OutcomeDS.StimType==40) = 50;
    OutcomeDS.StimVal(OutcomeDS.StimType==20) = 10;
    
    % Create a DS for previous 5 outcomes/selection
    Amt50_Outcome_idx = find(OutcomeDS.StimType==80 | OutcomeDS.StimType==40);
    Prev_5trls = [Amt50_Outcome_idx-1,Amt50_Outcome_idx-2,Amt50_Outcome_idx-3,Amt50_Outcome_idx-4,Amt50_Outcome_idx-5];
    Prev_5trls(Prev_5trls<1) = NaN;
    Prev_5trls_StimType = Prev_5trls;    
    Prev_5trls_StimVal = Prev_5trls;
    
    for iTrl = 1:size(Prev_5trls_StimVal,2);
        for iStim = 1:size(Prev_5trls_StimVal,1)
            if Prev_5trls(iStim,iTrl) > 0
                Prev_5trls_StimType(iStim,iTrl) = OutcomeDS.StimType(Prev_5trls(iStim,iTrl));
                Prev_5trls_StimVal(iStim,iTrl) = OutcomeDS.StimVal(Prev_5trls(iStim,iTrl));
            else
                continue
            end
        end
    end
    Prev_5trls_SrlNum_DS = mat2dataset(Prev_5trls,'VarNames',{'PrvTrl1','PrvTrl2','PrvTrl3','PrvTrl4','PrvTrl5'});    
    Prev_5trls_StimType_DS = mat2dataset(Prev_5trls_StimType,'VarNames',{'PrvTrl1','PrvTrl2','PrvTrl3','PrvTrl4','PrvTrl5'});
    Prev_5trls_StimVal_DS = mat2dataset(Prev_5trls_StimVal,'VarNames',{'PrvTrl1','PrvTrl2','PrvTrl3','PrvTrl4','PrvTrl5'});
    
    % Subject ID and Session
    TRBvarsDS.SubjID{iSubj} = SubjectIDs{iSubj};
    TRBvarsDS.SessRun{iSubj} = SessRun{iSubj};
    TRBvarsDS.Session{iSubj} = Session{iSubj};    
    
    % Number of Trials
    TRBvarsDS.N_Tot_Resp(iSubj) = size(OutcomeDS,1);
    TRBvarsDS.N50_nTrl(iSubj) = sum(OutcomeDS.StimType==80);
    TRBvarsDS.N10_nTrl(iSubj) = sum(OutcomeDS.StimType==60);
    TRBvarsDS.P50_nTrl(iSubj) = sum(OutcomeDS.StimType==40);
    TRBvarsDS.P10_nTrl(iSubj) = sum(OutcomeDS.StimType==20);
    
    % Extract RT and by selecting the first of the double button presses winin a single trial
    ChoiceDS = EventLogDS((EventLogDS.StimType==90),:);
    Ch_Events_N50 = ChoiceDS.EventNum(OutcomeDS.TrialNum(OutcomeDS.StimType==80));
    TRBvarsDS.N50_RT_M(iSubj) = mean(EventLogDS.TimeLag(Ch_Events_N50+1));
    Ch_Events_N10 = ChoiceDS.EventNum(OutcomeDS.TrialNum(OutcomeDS.StimType==60));  
    TRBvarsDS.N10_RT_M(iSubj) = mean(EventLogDS.TimeLag(Ch_Events_N10+1));
    Ch_Events_P50 = ChoiceDS.EventNum(OutcomeDS.TrialNum(OutcomeDS.StimType==40));  
    TRBvarsDS.P50_RT_M(iSubj) = mean(EventLogDS.TimeLag(Ch_Events_P50+1));
    Ch_Events_P10 = ChoiceDS.EventNum(OutcomeDS.TrialNum(OutcomeDS.StimType==20));  
    TRBvarsDS.P10_RT_M(iSubj) = mean(EventLogDS.TimeLag(Ch_Events_P10+1));
    
    TRBvarsDS.Loss_RT_M(iSubj) = mean([TRBvarsDS.N50_RT_M(iSubj);TRBvarsDS.N10_RT_M(iSubj)]);
    TRBvarsDS.Gain_RT_M(iSubj) = mean([TRBvarsDS.P50_RT_M(iSubj);TRBvarsDS.P10_RT_M(iSubj)]);
    TRBvarsDS.Amt50_RT_M(iSubj) = mean([TRBvarsDS.N50_RT_M(iSubj);TRBvarsDS.P50_RT_M(iSubj)]);
    TRBvarsDS.Amt10_RT_M(iSubj) = mean([TRBvarsDS.N10_RT_M(iSubj);TRBvarsDS.P10_RT_M(iSubj)]);
    
    % Selection Frequency for previous 'loss trials'
    TRBvarsDS.SF_50_LS_trials_1tr(iSubj) = sum(Prev_5trls_StimType_DS.PrvTrl1==80 | Prev_5trls_StimType_DS.PrvTrl1==60);
    SF_50_LS_trials_2tr_idx = intersect(find(Prev_5trls_StimType_DS.PrvTrl1==80 | Prev_5trls_StimType_DS.PrvTrl1==60), find(Prev_5trls_StimType_DS.PrvTrl2==80 | Prev_5trls_StimType_DS.PrvTrl2==60)); %#ok<*OR2>
    TRBvarsDS.SF_50_LS_trials_2tr(iSubj) = numel(SF_50_LS_trials_2tr_idx); 
    SF_50_LS_trials_3tr_idx = intersect(SF_50_LS_trials_2tr_idx, find(Prev_5trls_StimType_DS.PrvTrl3==80 | Prev_5trls_StimType_DS.PrvTrl3==60));
    TRBvarsDS.SF_50_LS_trials_3tr(iSubj) = numel(SF_50_LS_trials_3tr_idx); 
    SF_50_LS_trials_4tr_idx = intersect(SF_50_LS_trials_3tr_idx, find(Prev_5trls_StimType_DS.PrvTrl4==80 | Prev_5trls_StimType_DS.PrvTrl4==60));
    TRBvarsDS.SF_50_LS_trials_4tr(iSubj) = numel(SF_50_LS_trials_4tr_idx); 
    SF_50_LS_trials_5tr_idx = intersect(SF_50_LS_trials_4tr_idx, find(Prev_5trls_StimType_DS.PrvTrl5==80 | Prev_5trls_StimType_DS.PrvTrl5==60));
    TRBvarsDS.SF_50_LS_trials_5tr(iSubj) = numel(SF_50_LS_trials_5tr_idx); 
    
    % Selection Frequency for previous 'losing trend'
    TRBvarsDS.SF_50_LS_trend_2tr(iSubj) = sum((Prev_5trls_StimVal_DS.PrvTrl1 + Prev_5trls_StimVal_DS.PrvTrl2) < 0 );
    TRBvarsDS.SF_50_LS_trend_3tr(iSubj) = sum((Prev_5trls_StimVal_DS.PrvTrl1 + Prev_5trls_StimVal_DS.PrvTrl2 + Prev_5trls_StimVal_DS.PrvTrl3) < 0 );
    TRBvarsDS.SF_50_LS_trend_4tr(iSubj) = sum((Prev_5trls_StimVal_DS.PrvTrl1 + Prev_5trls_StimVal_DS.PrvTrl2 + Prev_5trls_StimVal_DS.PrvTrl3 + Prev_5trls_StimVal_DS.PrvTrl4) < 0 );    
    TRBvarsDS.SF_50_LS_trend_5tr(iSubj) = sum((Prev_5trls_StimVal_DS.PrvTrl1 + Prev_5trls_StimVal_DS.PrvTrl2 + Prev_5trls_StimVal_DS.PrvTrl3 + Prev_5trls_StimVal_DS.PrvTrl4 + Prev_5trls_StimVal_DS.PrvTrl5) < 0 );    
    
end

% Save the TRBvarsDS as a CSV file
%export(TRBvarsDS,'File','MGT_TRBvars_Ph4_n7603.csv','Delimiter',',');

%% END
