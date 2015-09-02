function [outdata]=baseline_erp_datastruct(indata,time_ms,base_window)

%baseline ERP data

%basewindow=[-100 0]

[~,n_chans,n_conds,n_subs]=size(indata);
outdata=zeros(size(indata));

[~,bl_s]=min(abs(time_ms-base_window(1)));
[~,bl_e]=min(abs(time_ms-base_window(2)));

for s=1:n_subs
    for cond=1:n_conds
        for chan=1:n_chans
            
            outdata(:,chan,cond,s)=indata(:,chan,cond,s) - ...
        mean(indata(bl_s:bl_e,chan,cond,s),1);
            
        end
    end
end

fprintf('Subject condition ERPs have been baselined from %d to %d ms\n', ...
    base_window(1), base_window(2))