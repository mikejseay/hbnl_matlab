function erpdata_out = filter_erpdata(erpdata_in, srate, cutoff)

%erpdata_in is timepts x channels x conditions x subjects

[EEG,n_conds,n_subs] = import_erpdata_eeglab_inline(erpdata_in,srate);

if length(cutoff) == 1
    EEG = pop_eegfiltnew(EEG, [], cutoff);
else
    EEG = pop_eegfiltnew(EEG, cutoff(1), cutoff(2) );
end

[n_channels,n_timepts,~]=size(EEG.data);

data = reshape(EEG.data,[n_channels,n_timepts,n_conds,n_subs]);
clear EEG

erpdata_out=permute(data,[2 1 3 4]);

end