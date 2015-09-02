function EEG = import_hdf_eeglab_inline(data_file,data,h1_struct)

%pull out the set name automatically, get the data, and sampling rate
[~,setname,~]=fileparts(data_file);
cnt_data=data';
cnt_srate=h1_struct.experiment_struct.rate; 

%import data
EEG = pop_importdata('setname',setname,'data',cnt_data,'dataformat','array',...
    'subject',setname(10:17),'session',setname(1:9),'nbchan',size(cnt_data,1),...
    'srate',cnt_srate);
EEG = eeg_checkset( EEG );

%add channel locations to remaining channels
EEG=pop_chanedit(EEG, 'lookup','/export/home/mike/matlab/eeglab12_0_2_6b/plugins/dipfit2.2/standard_BEM/elec/standard_1005.elc',...
    'load',{'/export/home/mike/matlab/coords/61chans_ns.ced' 'filetype' 'autodetect'},...
    'changefield',{47 'type' ''},'changefield',{47 'datachan' 1});
EEG = eeg_checkset( EEG );