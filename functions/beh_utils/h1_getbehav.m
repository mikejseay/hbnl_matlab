function [etable, h1_structout] = h1_getbehav(h1_file, trial_init_type)

% interpret h1_struct and param_struct to determine behavioral info,
% mainly response and reaction time

%trial_init_type=90;

if nargin<2
    trial_init_skip=true;
    trial_init_type=[];
else
    trial_init_skip=false;
end

%interpret h1_struct if it's just a file path
if ischar(h1_file)
    h1_file = read_hdf1_dataV7(h1_file);
end

stim_id=h1_file.case_struct.stim_id;
descriptor=h1_file.case_struct.descriptor;
response_id=h1_file.case_struct.response_id;

% fix problem with resp_code have values outside of [-1, 0, 1:8];
if any(~ismember(response_id, [-1, 0, 1:8]))
    response_id(~ismember(response_id, [-1, 0, 1:8])) = 0;
end

n_events=length(h1_file.event_struct.event_id);

%create a struct analogous to the CNT event log

event=1:n_events;
type=h1_file.event_struct.event_id;

trial_count=0;
trial=zeros(1,n_events);
type_descriptor=cell(1,n_events);
errant_response=false(1,n_events);
for ev=1:n_events
    if trial_init_skip
        trial_count=trial_count+1;
    elseif type(ev)==trial_init_type
        trial_count=trial_count+1;
    end
    trial(ev)=trial_count;
    if type(ev)==0  %some kind of response
        if ev>1     %this part catches errant responses
            if isempty(response_id(stim_id==type(ev-1))) || ...
                    response_id(stim_id==type(ev-1))==0
                type_descriptor{ev}='Errant Response';
                errant_response(ev)=true;
                continue
            end
        end
        type_descriptor{ev}='Response';
    else            %some kind of stimulus
        if ismember(type(ev),stim_id)
            type_descriptor{ev}=descriptor{stim_id==type(ev)};
        else
            type_descriptor{ev}='Unknown';
        end
    end
end

response_code=h1_file.event_struct.keypad_id;
event_time_s=h1_file.event_struct.event_time_offset;

%get the time between events
event_interval_ms=[0 round(diff(event_time_s)*1000)];

%get reaction times
rt=nan(1,n_events);
rt(response_code~=0 & ~errant_response) = ...
    event_interval_ms(response_code~=0 & ~errant_response);
%rts_correct=event_interval_ms(h1_struct.event_struct.event_id(1:end-1)==10 & ...
%    h1_struct.event_struct.keypad_id(2:end)==1);

cnth1_eventstruct=v2struct(event,trial,type,type_descriptor,...
    response_code,errant_response,event_time_s,event_interval_ms,rt);
cnth1_eventstruct=transposefields(cnth1_eventstruct);

etable=struct2table(cnth1_eventstruct);
h1_structout = h1_file;