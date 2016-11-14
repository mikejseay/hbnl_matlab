function bool = is_consecutive_boolean( ev_log, n )
% returns true if ev_log has at most n consecutive true or false values
% false otherwise

e = 1;
cc = 0;
max_e = length(ev_log);

while cc < n
    e = e + 1;
    if eq( ev_log(e-1), ev_log(e) )
        cc = cc + 1;
    else
        cc = 0;
    end
    if e == max_e
        bool = true;
        return
    end
end

bool = false;

end