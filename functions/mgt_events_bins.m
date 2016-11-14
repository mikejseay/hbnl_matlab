% script to create a random permutation of events for a behavioral task

n_events_perbin = 16;   % events per bin of trials
n_bins = 12;            % number of bins
n_consec_events = 5;    % number of consecutive events of one type

% proportion winning per bin uniformly random in [0.2 0.8]
%prop_win = make_binprobs(n_bins);
%prop_win = [0.5 .75 .25 .5 .25 .75 .5];
%prop_win = [0.5 .25 .75 .5 .75 .25 .5];
%prop_win = make_binprobs_frompool([0.2 0.5 0.8], n_bins);
prop_win = [.5 .8 .2 .5 .2 .8 .2 .5 .8 .5 .8 .2]';

%put space between bins?
space_bins = true;

%%

% this is the mapping between true and false and the stim codes
bool2stim_map = containers.Map([0 1],[1 4]);
bool2gosub_map = containers.Map([0 1],[1 2]);
%bool2gosub_map_pastcrit = containers.Map([0 1],[3 4]);
prob2cs_map = containers.Map([0.2 0.5 0.8],[10 20 30]);
prob2gosub_map = containers.Map([0.2 0.5 0.8],[200 100 200]);

%%

ev_log=[];

for bin=1:n_bins
    
    bin_log = make_evlog(n_events_perbin, prop_win(bin), n_consec_events);
    ev_log = [ev_log; bin_log];
    
end

%%

% convert this logical vector to an array formatted to be
% a stim sequence in the Stim2 Gentask program

cs_row      = {0, 'IMAGE', 800, 1000, 1506.67, 0, 0, 4, 90, 'c:\stimfiles\pfb\choice_stim.jpg'};
os_row      = {0, 'GOSUB', 999, 0, 0, 0, 0, 0, 0, 0};
check_row   = {0, 'IF', 'PC', 85, 100, 'GOTO', 0, 'NEXT', 0, 0};
blank_row   = {0, 'RESET', 0, 0, 0, 0, 0, 0, 0, 0};

n_events = length(ev_log);
n_columns = 10;
stim_array = cell(n_events, n_columns);
ed=0;
for ev=1:n_events
    
    %the choice stimulus
    ed=ed+1;
    stim_array(ed, :) = cs_row;
    
    %if first in block, override the label
    %if mod(ev, n_events_perbin) == 1
    %    stim_array{ed, 1} = 100 .* ( floor((ev-1)./n_events_perbin) + 1 );
    %end
    
    %override default correct response
    stim_array{ed, 8} = bool2stim_map(ev_log(ev));
    
    %override type code to indicate probability (?)
    % 0.2 = 10 (mostly left), 0.5 = 20 (equal), 0.8 = 30 (mostly right)
    blocknum = ceil( ev ./ n_events_perbin );
    stim_array{ed, 9} = prob2cs_map( prop_win( blocknum ) );
    
    %the outcome stimulus
    ed=ed+1;
    stim_array(ed, :) = os_row;
    
    %override default gosub code
    %these refer to subroutines with type codes to indicate probability
    % the gosubs are 0.5 =  101/102 (left/right), 0.2/0.8 = 201/202
    % (left/right). in turn, these have type codes of:
    % 0.5 = 40/50 (correct/incorrect), 0.2,0.8 = 60/70 (correct/incorrect)
    stim_array{ed, 3} = prob2gosub_map( prop_win( blocknum ) ) + ...
                            bool2gosub_map(ev_log(ev));
    
    blocknum = ceil( ev ./ n_events_perbin );
    stim_array{ed, 9} = prob2cs_map( prop_win( blocknum ) );
    
    % for trials 10-15, check to see if the % correct is > criterion, and if
    % so, GOTO the next block start
    if mod(ev, n_events_perbin) > 9
        ed=ed+1;
        stim_array(ed, :) = check_row;
        stim_array{ed, 7} = 100 .* ( floor((ev-1)./n_events_perbin) + 2 );
    end
    
    %space between bins
    if space_bins
        if mod(ev, n_events_perbin) == 0
            ed=ed+1;
            stim_array(ed, :) = blank_row;
            
            % override label
            stim_array{ed, 1} = 100 .* ( floor(ev./n_events_perbin) + 1 );
        end
    end
end

%% plot the probabilities of each block

evprob_vec = repmat(prop_win', [n_events_perbin 1]);

figure;
plot(evprob_vec(:),'b');
%hold on; plot(1 - evprob_vec(:),'m'); hold off;
axis([1 n_events 0 1]);
lh(1) = xlabel('Trial #');
lh(2) = ylabel('Probability of Right-Hand Deck Winning');
set(lh, 'fontsize', 24);
set(gca, 'fontsize', 24);