function ev_out = fix_consecutive_boolean( ev_in, n )
% create a version of a boolean vector with same number of true and false
% elements but with a maximum of n consecutive elements of each type

max_e = length(ev_in);  % total events
n_true = sum(ev_in);    % number of true events at input

ev_tmp = ev_in;

while ~is_consecutive_boolean(ev_tmp, n) %big while
    
e = 1;  % event counter
cc = 0; % consecutivity counter

while cc < n
    e = e + 1;
    if eq( ev_tmp(e-1), ev_tmp(e) )
        cc = cc + 1;
    else
        cc = 0;
    end
end

% pick a random different opposite logical event to invert
other_ev = randperm(max_e, 1);
% repick if the same event index or not opposite
while other_ev == e || eq( ev_tmp(e),  ev_tmp(other_ev) )
    other_ev = randperm(max_e, 1);
end

% invert them both
ev_tmp(e) = ~ev_tmp(e);
ev_tmp(other_ev) = ~ev_tmp(other_ev);

% verify number of truthes (should not happen but just in case)
if sum(ev_tmp) < n_true
    maketrue_ev = randperm(max_e, 1);
    while ev_tmp(maketrue_ev) % if true repick
        maketrue_ev = randperm(max_e, 1);
    end
    ev_tmp(maketrue_ev) = true;
elseif sum(ev_tmp) > n_true
    makefalse_ev = randperm(max_e, 1);
    while ~ev_tmp(makefalse_ev) % if false repick
        makefalse_ev = randperm(max_e, 1);
    end
    ev_tmp(makefalse_ev) = false;
end

end % end of big while

ev_out = ev_tmp;

end