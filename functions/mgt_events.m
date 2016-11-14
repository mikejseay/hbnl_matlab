% script to create a random permutation of events for a behavioral task

n_events = 192;         % total events
prop_win = 0.75;        % proportion "true"
n_consec_events = 5;    % number of consecutive events of one type

% this is the mapping between true and false and the stim codes
bool2code_map = containers.Map([true false],[999 5999]);

%%

% create logical vector of events (1 = win, 0 = loss)

n_wins = floor(n_events * prop_win);
n_losses = n_events - n_wins;

win_inds = sort(randperm(n_events, n_wins));

ev_log = false(n_events, 1);
ev_log(win_inds) = true;

ev_log2 = fix_consecutive_boolean(ev_log, n_consec_events);

%%

% convert this logical vector to an array formatted to be
% a stim sequence in the Stim2 Gentask program

stimcode_array = zeros(n_events*2, 1);
ed=0;
for ev=1:n_events
    ed=ed+1;
    stimcode_array(ed) = 800;
    ed=ed+1;
    stimcode_array(ed) = bool2code_map(ev_log2(ev));
end