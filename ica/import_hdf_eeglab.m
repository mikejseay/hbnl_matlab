function EEG = import_hdf_eeglab(data_file,chans_to_remove)

if nargin<2
    %these are the defaultly "bad" chans
    %32=HEOG,63="blank",64=VEOG
    chans_to_remove=[32 63 64];
end

%read the HDF file
h1_struct = read_hdf1_dataV7(data_file);

%make certain vars global
%global cnt_data 

%pull out the set name automatically, get the data, and sampling rate
[pathstr,setname,~]=fileparts(data_file);
cnt_data=h1_struct.data_struct.hdf1_cnt_data(:,1:end-19);
cnt_srate=h1_struct.experiment_struct.rate; 

%start EEGLAB
%[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%import data
EEG = pop_importdata('setname',setname,'data',cnt_data,'dataformat','array',...
    'subject',setname(10:17),'session',setname(1:9),'nbchan',size(cnt_data,1),...
    'srate',cnt_srate);

%pull data into eeglab and check it
%eeglab redraw

%remove extra 3 channels
EEG = pop_select( EEG,'nochannel',chans_to_remove);
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );

%add channel locations to remaining channels
EEG=pop_chanedit(EEG, 'lookup','/export/home/mike/matlab/eeglab13_4_4b/plugins/dipfit2.3/standard_BEM/elec/standard_1005.elc',...
    'load',{'/export/home/mike/matlab/origin/coords/61chans_ns.ced' 'filetype' 'autodetect'},...
    'changefield',{47 'type' ''},'changefield',{47 'datachan' 1});
%EEG = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );


%create an event array
type=h1_struct.trial_struct.case_num;
latency=h1_struct.trial_struct.time_offset;
rt=h1_struct.trial_struct.response_latency;
resp=h1_struct.trial_struct.response_id;
trial_events=[type;latency;rt;resp]';

%import events
EEG = pop_importevent(EEG, 'event', trial_events, 'fields',...
    {'type', 'latency', 'rt', 'resp'}, 'timeunit', 1);

%redraw eeglab
eeglab redraw