function [erpdata,peakmat] = coh_pickERPpeaks (imp, scl, s_inds_g, erpdata, opt_logic, lp_cutoff, base_window)
% filter, baseline, and pick peaks in erpdata

if nargin < 7
    base_window=[-100 0];
end
if nargin < 6
    lp_cutoff = 16;
end
if nargin < 5
    opt_logic = [true true];
end

% filter
if opt_logic(1)
    erpdata=filter_erpdata(erpdata,imp.erptimerate,lp_cutoff);
end

% baseline
if opt_logic(2)
    erpdata=baseline_erp_datastruct(erpdata,scl.t_ms_erp, ...
        base_window);
end

% pick peaks

%specify chan
peakchan=7;

%specify window in which peak may be found
pkwin_s_ms=250;
pkwin_e_ms=550;

%pause for quality control?
dopause=false;

% P3
[~,pkwin_s]=min(abs(scl.t_ms_erp-pkwin_s_ms));
[~,pkwin_e]=min(abs(scl.t_ms_erp-pkwin_e_ms));
peakmat=zeros(size(s_inds_g,1),imp.maxconds,2);
%f=figure;
for s=1:size(s_inds_g,1)
if s_inds_g(s,1)
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
    peakmat(s,cond,:)=[val,scl.t_ms_erp(pkwin_s+ind-1)];
end
if dopause
    pause
end
%clf
end
end
%close(f)

end