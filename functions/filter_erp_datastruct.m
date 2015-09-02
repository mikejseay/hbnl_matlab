function [outdata]=filter_erp_datastruct(indata,srate,lowpass_band)

%lowpass filter ERP data

%srate=256;
%lowpass_band=16;

[~,n_chans,n_conds,n_subs]=size(indata);
outdata=zeros(size(indata));

for s=1:n_subs
    for cond=1:n_conds
        for chan=1:n_chans
            
            outdata(:,chan,cond,s)=eegfilt(indata(:,chan,cond,s)',srate,0,lowpass_band);
            
        end
    end
end

fprintf('Subject condition ERPs have been lowpass filtered at %d Hz.\n',lowpass_band)