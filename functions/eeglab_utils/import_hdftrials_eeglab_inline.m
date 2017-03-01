function EEG = import_hdftrials_eeglab_inline(data_file,data,h1_struct_orrate)

if isstruct(h1_struct_orrate)
    cnt_srate=h1_struct_orrate.experiment_struct.rate; 
elseif isnumeric(h1_struct_orrate)
    cnt_srate=h1_struct_orrate;
end

%pull out the set name automatically, get the data, and sampling rate
[~,setname,~]=fileparts(data_file);
cnt_data=permute(data,[3 1 2]);

%import data, must be channels x timepts x trials
EEG = pop_importdata('setname',setname,'data',cnt_data,'dataformat','array',...
    'subject',setname(10:17),'session',setname(1:9),'nbchan',size(cnt_data,1),...
    'srate',cnt_srate);
EEG = eeg_checkset( EEG );

%add channel locations to remaining channels
EEG=pop_chanedit(EEG, 'lookup','/active_projects/matlab_common/standard_1005.elc',...
    'load',{'/active_projects/matlab_common/61chans_ns.ced' 'filetype' 'autodetect'},...
    'changefield',{47 'type' ''},'changefield',{47 'datachan' 1});
EEG = eeg_checkset( EEG );