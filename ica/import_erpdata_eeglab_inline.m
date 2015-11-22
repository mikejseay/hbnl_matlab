function [EEG,n_conds,n_subs] = import_erpdata_eeglab_inline(data,cnt_srate)

%at this point, the data is timepts x channels x conditions x subjects
[n_timepts,n_channels,n_conds,n_subs]=size(data);

data2=reshape(data,[n_timepts n_channels n_conds*n_subs]);
clear data
%now the data is timepts x channnels x "trials"

%change to channels x timepts x trials
cnt_data=permute(data2,[2 1 3]);
clear data2

%import data, must be channels x timepts x trials
EEG = pop_importdata('setname',' ','data',cnt_data,'dataformat','array',...
    'subject',' ','session',' ','nbchan',size(cnt_data,1),...
    'srate',cnt_srate);
EEG = eeg_checkset( EEG );

%add channel locations to remaining channels
%EEG=pop_chanedit(EEG, 'lookup','A:/matlab/eeglab12_0_2_6b/plugins/dipfit2.2/standard_BEM/elec/standard_1005.elc',...
%    'load',{'A:/matlab/origin/coords/61chans_ns.ced' 'filetype' 'autodetect'},...
%    'changefield',{47 'type' ''},'changefield',{47 'datachan' 1});
%EEG = eeg_checkset( EEG );