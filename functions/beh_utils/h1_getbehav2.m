function etable = h1_getbehav2(h1_file, expstruct)
% interpret cnt event log given expstruct, a structure containing
% experimental attributes

%interpret h1_struct if it's just a file path
if ischar(h1_file)
    h1_file = read_hdf1_dataV7(h1_file);
end

% rename attributes of the expstruct that describes the experiment

stim_id         = expstruct.ttl_codes;
descriptor      = expstruct.ttl_descriptors;
response_id     = expstruct.correct_resps;
response_win    = expstruct.resp_windows;
min_resptime    = expstruct.min_resptime;
trial_init_type = expstruct.init_trial;

%create a struct analogous to the CNT event log

n_events        = length(h1_file.event_struct.event_id);

event           = 1:n_events;
type            = h1_file.event_struct.event_id;
resp_code       = h1_file.event_struct.keypad_id;
event_time_s    = h1_file.event_struct.event_time_offset;

% fix problem with resp_code have values outside of [-1, 0, 1:8];
if any(~ismember(resp_code, [-1, 0, 1:8]))
    resp_code(~ismember(resp_code, [-1, 0, 1:8])) = 0;
end

trial           = zeros(1,n_events);
type_descriptor = cell(1,n_events);

% only apply to responses
[ errant, early, late, correct ] = deal( false(1,n_events) );

% step through event log
trial_count = 0;
for ev=1:n_events
    if ismember(type(ev), trial_init_type)
        trial_count=trial_count+1;
    end
    trial(ev) = trial_count;
    if type(ev)==0  %some kind of response
    if ev>1
       
        % this part catches errant responses, defined as responses that
        % were made at the wrong part of the experiment (when no
        % response was called for) or immediately following a previous
        % response (double press). these responses are not considered
        % for the purpose of "hit rates" or "false alarms" (i.e.
        % accuracy).
        tmp_prevtype = stim_id == type(ev-1);
        if isempty( response_id(tmp_prevtype) ) || ...
                response_id(tmp_prevtype) == -1
            type_descriptor{ev} = 'rsp_err';
            errant(ev) = true;
            continue % is not considered early/late or correct/incorrect
        end
        
        % this next part catches early / late responses. early responses
        % are usually considered as "no responses," i.e. the subject did
        % not actually react so the event should not be interpreted. late
        % responses may be considered errors by some. this algorithm allows
        % users to interpret the event table as they would like but still
        % goes on to interpret its correctness
        tmp_rt = (event_time_s(ev) - event_time_s(ev-1))*1000;
        if tmp_rt > response_win(tmp_prevtype)
            late(ev) = true;
            early(ev) = false;
        elseif tmp_rt < min_resptime
            early(ev) = true;
            late(ev) = false;
            type_descriptor{ev}='rsp_early';
            continue
        end
        
        % this next part interprets the response as correct / incorrect
        % note that it could have been late
        if resp_code(ev) == response_id(tmp_prevtype)
            type_descriptor{ev}='rsp_correct';
            correct(ev) = true;
            continue
        else
            type_descriptor{ev}='rsp_incorrect';
            correct(ev) = false;
            continue
        end
        
    end
    type_descriptor{ev}='rsp_unknown';
    
    else            %some kind of stimulus
        
        % if a code for this task
        if ismember(type(ev), stim_id)
            tmp_currtype = stim_id==type(ev);
            type_descriptor{ev} = descriptor{tmp_currtype};
            
            % if the last event
            if ev+1 > n_events
                % if no response was required, it's correct
                if ismember(response_id(tmp_currtype), [0 -1])
                    correct(ev) = true;
                end
            else % if not the last event
                % if the following response was correct
                if resp_code(ev+1) == response_id(tmp_currtype)
                    correct(ev) = true;
                else %otherwise, it's not correct
                    correct(ev) = false;
                end
            end
        else % if not a code for this task
            type_descriptor{ev} = 'stm_unknown';
        end
    end
end

%get the time between events
event_interval_ms = [ 0 diff(event_time_s)*1000 ];

%get reaction times
rt = nan(1,n_events);
rt(resp_code~=0) = event_interval_ms(resp_code~=0);
rt(errant==true) = NaN;

h1_eventstruct = v2struct(event, trial, type, type_descriptor,...
    correct, rt, resp_code, errant, early, late, event_time_s, ...
    event_interval_ms);
h1_eventstruct = transposefields(h1_eventstruct);

etable = struct2table(h1_eventstruct);