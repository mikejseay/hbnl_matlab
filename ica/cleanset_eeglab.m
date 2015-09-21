function EEGOUT = cleanset_eeglab(EEGIN)

arg_flatline=[]; %'off';
arg_highpass=[];
arg_channel=[]; %'off';
arg_noisy=[]; %'off';
arg_burst=10;
arg_window='off';
                            
EEGOUT=clean_rawdata(EEGIN,arg_flatline,arg_highpass,...
    arg_channel,arg_noisy,arg_burst,arg_window);