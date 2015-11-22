function erpdata_out = filter_erpdata(erpdata_in,srate,lp_cutoff)

[EEG,n_conds,n_subs] = import_erpdata_eeglab_inline(erpdata_in,srate);
EEG = pop_eegfiltnew(EEG, [], lp_cutoff);

[n_channels,n_timepts,~]=size(EEG.data);

data = reshape(EEG.data,[n_channels,n_timepts,n_conds,n_subs]);
clear EEG

erpdata_out=permute(data,[2 1 3 4]);

end