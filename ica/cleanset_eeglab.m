function EEGOUT = cleanset_eeglab(EEGIN,hpfilt_transition_bands)

if nargin<2
    arg_highpass=[];
else
    arg_highpass=hpfilt_transition_bands;
end

arg_flatline=[]; %'off';
arg_channel=[]; %'off';
arg_noisy=[]; %'off';
arg_burst=[];
arg_window='off'; %off because it messes up the trial indexing scheme
 %as a replacement for this step, we reject trials based on high voltage
                            
EEGOUT=clean_rawdata(EEGIN,arg_flatline,arg_highpass,...
    arg_channel,arg_noisy,arg_burst,arg_window);