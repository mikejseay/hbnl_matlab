function chan_locs=pass_chan_locs(chan_thing)

if isscalar(chan_thing)
    %if coords spec is a number, load the file with that number of channels
    switch chan_thing
        case 31
            load('/export/home/mike/matlab/origin/coords/31chans_ns.mat')
        case 61
            load('/export/home/mike/matlab/origin/coords/61chans_ns.mat')
    end
elseif ischar(chan_thing)
    %if a string, load that location file
    load(chan_thing)
elseif isstruct(chan_thing)
    %if a structure, set that as the chan_locs
    chan_locs=chan_thing;
else
    error('Channel locations incorrectly specified.')
end

end