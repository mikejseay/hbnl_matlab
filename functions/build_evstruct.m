function [ev, urev] = build_evstruct( n_samps, srate, seg_s, name )

seglength = round(srate * seg_s);

lats = 1:seglength:n_samps;
n_events = length(lats);

type = cell(n_events, 1);
[type{:}]=deal(name);

urev = struct('type', type, 'latency', num2cell(lats)');

urev_inds = 1:n_events;
ev = struct('type', type, 'latency', num2cell(lats)', 'urevent', num2cell(urev_inds)');

end