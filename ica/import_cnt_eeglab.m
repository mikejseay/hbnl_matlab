%define the cnt file, dat file, and generic epoch length\
%data_file='/vol01/active_projects/crystalann/washu/ns374/50723006/ant_6_d1_50723006_32.cnt';
data_file='/raw_data/neuroscan/suny/ns164/a0001071/vp3_4_a1_a0001071.cnt';
%data_file='/vol01/raw_data/staging/iowa/ns316/30172356/ant_6_d1_30172356_32.cnt';
%'C:\Users\Mike\Desktop\ica test\ern_9_a1_40001009_32.cnt';
%dat_file='C:\Users\Mike\Desktop\ica test\ern_9_a1_40001009.dat';
%epoch_length=1506;

%start EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%import data
EEG = pop_loadcnt(data_file, 'dataformat', 'auto', 'memmapfile', '');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
EEG = eeg_checkset( EEG );

%remove extra 3 channels
%EEG = pop_select( EEG,'nochannel',[32 63 64] );
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
%EEG = eeg_checkset( EEG );

%add channel locations to remaining channels
%EEG=pop_chanedit(EEG, 'lookup','A:\\matlab\\eeglab13_4_4b\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc','load',...
%    {'A:\\matlab\\coords\\61chans_ns.ced' 'filetype' 'autodetect'},'changefield',{47 'type' ''},'changefield',{47 'datachan' 1});
%[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%filter
EEG = pop_eegfiltnew(EEG, 1, 50, 1650, 0, [], 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );

%reference to average
%EEG = pop_reref( EEG, []);
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 

%down-sample timepoints


%clean data using ASR?
%EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5);
%EEG = eeg_checkset( EEG );

%extract epochs
%EEG = pop_epoch( EEG, {  }, [-0.5 1], 'newname', 'CNT file epochs', 'epochinfo', 'yes');
%[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
%EEG = eeg_checkset( EEG );

%import behavioral data
%EEG = pop_loaddat( EEG, dat_file, epoch_length);
%[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%redraw eeglab
eeglab redraw