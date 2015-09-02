% pick peaks in erpdata

%% filter (run only once)

filter_lowpass=16;
erpdata=filter_erp_datastruct(erpdata,param_struct.rate, ...
    filter_lowpass);

%% baseline (run only once)

baseline_window_ms=[-100 0];
erpdata=baseline_erp_datastruct(erpdata,scl.t_ms, ...
    baseline_window_ms);

%% look at peaks

%specify chan
peakchan=7;

%specify window in which peak may be found
pkwin_s_ms=250;
pkwin_e_ms=550;

%pause for quality control?
%dopause=false;
dopause=false;

%
[~,pkwin_s]=min(abs(scl.t_ms-pkwin_s_ms));
[~,pkwin_e]=min(abs(scl.t_ms-pkwin_e_ms));
peakmat=zeros(size(s_inds_g,1),imp.maxconds,2);
f=figure;
for s=1:size(s_inds_g,1)
if s_inds_g(s,9)
subplot_dummy=0;
for cond=1:imp.maxconds
    erp_plot_data=erpdata(:,peakchan,cond,s);
    [ind,val]=peakfinder(erp_plot_data(pkwin_s:pkwin_e),[],[],[],false);
    if isempty(val)
        [val,ind]=max(erp_plot_data(pkwin_s:pkwin_e));
    end
    if length(val) > 1 || length(ind) > 1
        [val,tempind]=max(val);
        ind=ind(tempind);
        clear tempind
    end
    %
    if false
    subplot_dummy=subplot_dummy+1;
    subplot(sp_d(1),sp_d(2),subplot_dummy)
    plot(erp_plot_data); hold on;
    plot(pkwin_s+ind-1,val,'rx');
    axis auto; set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); grid on;
    vline(scl.t_zero,'k--'); vline(pkwin_s,'r'); vline(pkwin_e,'r');
    title(sprintf('S%d',s)); hold off;
    end
    %
    peakmat(s,cond,:)=[val,pkwin_s+ind-1];
end
if dopause
    pause
end
clf
end
end
close(f)

clear bl_s_ms bl_e_ms bl_s bl_e chan cond s peakchan pkwin_s_ms pkwin_e_ms ...
    dopause pkwin_s pkwin_e f erp_plot_data ind val subplot_dummy ...
    baseline_window_ms filter_lowpass