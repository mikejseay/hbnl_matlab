function EEGOUT = cleanset_eeglab(EEGIN, arg_highpass, arg_burst)

if nargin<3
    arg_burst = 'off';
end

if nargin<2
    arg_highpass = [];
end

arg_flatline = [];
arg_channel = [];
arg_noisy = [];
arg_window = 'off'; %off because it messes up the trial indexing scheme
% (it rejects windows of data)
                            
EEGOUT=clean_rawdata(EEGIN,arg_flatline,arg_highpass,...
    arg_channel,arg_noisy,arg_burst,arg_window);