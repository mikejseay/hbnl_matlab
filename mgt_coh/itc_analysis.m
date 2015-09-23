%% check ERP peak data

figure; subplot_dummy=0;
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    sp(subplot_dummy)=subplot(sp_d(1),sp_d(2),subplot_dummy);
    for group=pp.chosen_g
        scatter(peakmat(s_inds_g(:,group),cond,2),...
            peakmat(s_inds_g(:,group),cond,1),scl.g_color{group}); hold on;
    end
    hold off; set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms);
    axis tight; title(sprintf('%s',scl.cond_label{cond}));
    xlabel('Peak Latency (ms)'); ylabel('Peak Amplitude (uV)');
end
linkaxes(sp)
clear_plotassistvars

%% plot ERPS with groups superimposed

h_line=zeros(length(pp.plotn_cond),length(pp.chosen_g));
for chan=pp.chosen_chan(pp.plotn_chan)
    figure;
    subplot_dummy=0;
    for cond=pp.plotn_cond
        subplot_dummy=subplot_dummy+1;
        sp(cond)=subplot(sp_d(1),sp_d(2),subplot_dummy);
        for group=pp.chosen_g(pp.plotn_g)
            if cond==imp.maxconds+1
                erp_plot_data = mean(mean(erpdata(:,chan,pp.cond_diff{1},s_inds_g(:,group)),3),4) - ...
            mean(mean(erpdata(:,chan,pp.cond_diff{2},s_inds_g(:,group)),3),4);
                erp_plot_data_std = std(mean(erpdata(:,chan,pp.cond_diff{1},s_inds_g(:,group)),3) - ...
            mean(erpdata(:,chan,pp.cond_diff{2},s_inds_g(:,group)),3),0,4)/sqrt(sum(s_inds_g(:,group)));
            else
                erp_plot_data=mean(erpdata(:,chan,cond,s_inds_g(:,group)),4);
                erp_plot_data_std=std(erpdata(:,chan,cond,s_inds_g(:,group)),0,4)/sqrt(sum(s_inds_g(:,group)));
            end
            %plot(ero_plot_data,scl.g_color{group}); hold on; %,'Color',scl.s_color(pp.chosen_s,:)); hold on
            h=shadedErrorBar(1:imp.maxtimepts,erp_plot_data,erp_plot_data_std,scl.g_color{group}); hold on;
            h_line(cond,group)=h.mainLine;
        end
        axis([scl.t_start scl.t_end -6 10]);
        vline(scl.t_zero,'k--'); hold off;
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
        title([scl.chan_label{chan},'/',scl.cond_label{cond}])
    end
tightfig;
linkaxes(sp(1:end-1))
end
clear_plotassistvars

%% plot ERPs with conditions superimposed

cond_colors={'r','g','b','m','r--','g--','b--','m--','k--','c'};
for chan=pp.chosen_chan(pp.plotn_chan)
    figure;
    group_dummy=0;
    for group=pp.chosen_g(pp.plotn_g)
        group_dummy=group_dummy+1;
        sp(group_dummy)=subplot(1,length(pp.plotn_g),group_dummy);
        for cond=pp.plotn_cond
            if cond==imp.maxconds+1
                erp_plot_data = mean(mean(erpdata(:,chan,pp.cond_diff{1},s_inds_g(:,group)),3),4) - ...
            mean(mean(erpdata(:,chan,pp.cond_diff{2},s_inds_g(:,group)),3),4);
                erp_plot_data_std = std(mean(erpdata(:,chan,pp.cond_diff{1},s_inds_g(:,group)),3) - ...
            mean(erpdata(:,chan,pp.cond_diff{2},s_inds_g(:,group)),3),0,4)/sqrt(sum(s_inds_g(:,group)));
            else
                erp_plot_data=mean(erpdata(:,chan,cond,s_inds_g(:,group)),4);
                erp_plot_data_std=std(erpdata(:,chan,cond,s_inds_g(:,group)),0,4)/sqrt(sum(s_inds_g(:,group)));
            end
            plot(erp_plot_data,cond_colors{cond}); hold on;
        end
        axis tight
        vline(scl.t_zero,'k--'); hold off;
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
        title([scl.chan_label{chan},'/',scl.g_label{group}])
        clickableLegend(scl.cond_label(pp.plotn_cond))
    end
tightfig;
linkaxes(sp)
end
clear_plotassistvars

%% plot topography of ERPs

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

erp_topo_scale=[-5 3];
erp_diff_limits=[-2 2];
for group=pp.chosen_g(pp.plotn_g)
figure; subplot_dummy=0;
overtitle=sprintf('Topography of ERPs / %s',scl.g_label{group});
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)
    [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    sp(cond)=subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy);
    if cond==imp.maxconds+1
        erp_topo_data = mean(mean(mean(erpdata(t_start:t_end,:,pp.cond_diff{1},s_inds_g(:,group)),1),3),4) - ...
    mean(mean(mean(erpdata(t_start:t_end,:,pp.cond_diff{2},s_inds_g(:,group)),1),3),4);
        topoplot(erp_topo_data,chan_locs,'maplimits',[erp_diff_limits(1) erp_diff_limits(2)],'electrodes','off','colormap',pp.cmap);
    else
        erp_topo_data=mean(mean(erpdata(t_start:t_end,:,cond,s_inds_g(:,group)),1),4);
        topoplot(erp_topo_data,chan_locs,'maplimits',[erp_topo_scale(1) erp_topo_scale(2)],'electrodes','off','colormap',pp.cmap);
    end
    %title(sprintf('%s, %d - %d ms, %s',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),scl.g_label{group}))
    end
    colorbar;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
tightfig;
end
clear_plotassistvars

%% look at EROs in individual bands / superimpose groups

%h_line=zeros(length(pp.plotn_cond),length(pp.chosen_g));
for chan=pp.chosen_chan(pp.plotn_chan)
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    figure;
    subplot_dummy=0;
    for cond=pp.plotn_cond
        subplot_dummy=subplot_dummy+1;
        sp(cond)=subplot(sp_d(1),sp_d(2),subplot_dummy);
        for group=pp.chosen_g(pp.plotn_g)
            if cond==imp.maxconds+1
                ero_plot_data = mean(mean(mean(wave_totdata(:,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),3),4),5) - ...
            mean(mean(mean(wave_totdata(:,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),3),4),5);
                ero_plot_data_std = std(mean(mean(wave_totdata(:,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),3),4) - ...
            mean(mean(wave_totdata(:,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),3),4),0,5)/sqrt(sum(s_inds_g(:,group)));
            else
                ero_plot_data=mean(mean(wave_totdata(:,chan,f_end:f_start,cond,s_inds_g(:,group)),3),5);
                ero_plot_data_std=std(mean(wave_totdata(:,chan,f_end:f_start,cond,s_inds_g(:,group)),3),0,5)/sqrt(sum(s_inds_g(:,group)));
            end
            %plot(ero_plot_data,scl.g_color{group}); hold on; %,'Color',scl.s_color(pp.chosen_s,:)); hold on
            h=shadedErrorBar(1:imp.maxtimepts,ero_plot_data,ero_plot_data_std, scl.g_color{group}); hold on;
            %h_line(cond,group)=h.mainLine;
        end
        vline(scl.t_zero,'k--'); hold off;
        axis tight
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
        title([scl.chan_label{chan},'/',scl.cond_label{cond},'/',num2str(pp.f_start_hz(freq_range)),'-',num2str(pp.f_end_hz(freq_range)),' Hz'])
    end
tightfig;
linkaxes(sp(1:end-1))
end
end
clear_plotassistvars

%% look at EROS in individual bands / superimpose conditions

ero_yax=[-20 80];

cond_colors={'r','g','b','m','r--','g--','b--','m--','k--','c'};
for chan=pp.chosen_chan(pp.plotn_chan)
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    figure;
    overtitle=sprintf('%1.1f - %1.1f at %s',pp.f_start_hz(freq_range), ...
        pp.f_end_hz(freq_range),scl.chan_label{chan});
    group_dummy=0;
    for group=pp.chosen_g(pp.plotn_g)
        group_dummy=group_dummy+1;
        sp(group_dummy)=subplot(1,length(pp.plotn_g),group_dummy);
        for cond=pp.plotn_cond
            if cond==imp.maxconds+1
                ero_plot_data = mean(mean(mean(wave_totdata(:,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),3),4),5) - ...
            mean(mean(mean(wave_totdata(:,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),3),4),5);
                ero_plot_data_std = std(mean(mean(wave_totdata(:,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),3),4) - ...
            mean(mean(wave_totdata(:,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),3),4),0,5)/sqrt(sum(s_inds_g(:,group)));
            else
                ero_plot_data=mean(mean(wave_totdata(:,chan,f_end:f_start,cond,s_inds_g(:,group)),3),5);
                ero_plot_data_std=std(mean(wave_totdata(:,chan,f_end:f_start,cond,s_inds_g(:,group)),3),0,5)/sqrt(sum(s_inds_g(:,group)));
            end
            plot(ero_plot_data,cond_colors{cond}); hold on;
        end
        axis([1 scl.t_end ero_yax(1) ero_yax(2)]);
        vline(scl.t_zero,'k--'); hold off;
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
        title(scl.g_label{group})
        clickableLegend(scl.cond_label(pp.plotn_cond),'Location','best')
    end
tightfig;
plottitle(overtitle);
%linkaxes(sp);
end
end
clear_plotassistvars

%% image ERO (ERSP) in time-freq at a chosen channel

pp.figdum=pp.figdum_init;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
for group=pp.chosen_g(pp.plotn_g)
for chan=pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(sp_d(1),sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        ero_plot_data=squeeze(mean(mean(wave_totdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),4),5)-...
            mean(mean(wave_totdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),4),5));
        %ero_plot_data=squeeze(mean(wavelet_tot(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),5)-...
        %    mean(wavelet_tot(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),5));
    else
        %ero_plot_data=squeeze(mean(wave_totdata(:,chan,:,cond,s_inds_g(:,group)),5));
        ero_plot_data=squeeze(mean(wave_totdata(:,chan,:,cond,s_inds_g(:,group)),5)) - ... %this one is subtractively "baseline-normalized"
            repmat(squeeze(mean(mean(wave_totdata(1:scl.t_zero,chan,:,cond,s_inds_g(:,group)),1),5)),1,imp.maxtimepts)';
        %ero_plot_data=squeeze(mean(wave_totdata(:,chan,:,cond,s_inds_g(:,group)),5)) ./ ... %this one is divisively "baseline-normalized"
        %    repmat(squeeze(mean(mean(wave_totdata(1:scl.t_zero,chan,:,cond,s_inds_g(:,group)),1),5)),1,imp.maxtimepts)';
    end
    contourf(fliplr(ero_plot_data)',pp.n_contour)
    shading flat; colormap(pp.cmap);
    %imagesc(ero_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,cond,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); ylabel('Frequency (Hz)');
    grid on;
    title(['ERO at ',scl.chan_label{chan},' : ',scl.g_label{group},' : ',scl.cond_label{cond}])
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
end
end
%
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
for splot=1:subplot_dummy
    if splot==subplot_dummy
        subplot(sp_d(1),sp_d(2),splot); caxis([c_diff(1) c_diff(2)]);
    else
        subplot(sp_d(1),sp_d(2),splot); caxis([c(1) c(2)]);
    end
    colorbar;
end
tightfig;
end
%distFig('s','ext','transpose',true);
clear_plotassistvars

%% plot ERO as a topographic plot over the head

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

ero_topo_scale=[0 15];
ero_diff_limits=[-8 2];

v=zeros(length(pp.plotn_cond),length(pp.chosen_g),length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms),2);
for group=pp.chosen_g(pp.plotn_g)
for freq_range=pp.plotn_f
figure; subplot_dummy=0;
overtitle=sprintf('Topography of EROs in %1.1f - %1.1f Hz / %s',pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.g_label{group});
%convert scl.freqs to pts
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)
        %convert times to points
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        %plot
        subplot_dummy=subplot_dummy+1;
        subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy)
        if cond==imp.maxconds+1
            topo_data=squeeze(mean(mean(mean(mean(wave_totdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4),5)) - ...
                squeeze(mean(mean(mean(mean(wave_totdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4),5));
            topoplot(topo_data,chan_locs,'maplimits',[ero_diff_limits(1) ero_diff_limits(2)],'electrodes','off','colormap',pp.cmap);
        else
            %topo_data=squeeze(mean(mean(mean(wave_totdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5));
            topo_data=squeeze(mean(mean(mean(wave_totdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5)) -...
                squeeze(mean(mean(mean(wave_totdata(1:scl.t_zero,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5));
            topoplot(topo_data,chan_locs,'maplimits',[ero_topo_scale(1) ero_topo_scale(2)],'electrodes','off','colormap',pp.cmap);
        end
        v(cond,group,freq_range,win,1) = min(topo_data); v(cond,group,freq_range,win,2) = max(topo_data);
        %title(sprintf('%s, %d - %d ms, %1.1f - %1.1f Hz, %s',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.g_label{group}))
    end
    colorbar;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
%tightfig;
end
end
c(1)=min(min(min(min(v(1:end-1,:,:,:,1))))); c(2)=max(max(max(max(v(1:end-1,:,:,:,2)))));
c_diff(1)=min(min(min(v(end,:,:,:,1)))); c_diff(2)=max(max(max(v(end,:,:,:,2))));
clear_plotassistvars

%% scatter ERO with RT measures

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='RT Variance (ms)'; %Median, variance, std also interesting
y_plotlabel='ERO Power';

behav_data=varRT;

ero_axes=[min(min(behav_data)) 60000 -20 20]; %max(max(behav_data))
ero_axes_diff=[min(min(behav_data)) max(max(behav_data)) -2 5];

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_start_hz(pp.plotn_f)),pp.maxwin);
for chan=pp.chosen_chan(pp.plotn_chan)
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    figure; subplot_dummy=0;
    overtitle=sprintf('ERO Power over age at %s in %1.1f - %1.1f Hz', ...
        scl.chan_label{chan},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
    for cond=pp.plotn_cond
    for win=1:pp.maxwin
    [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
    %    
    for group=pp.chosen_g(pp.plotn_g)
        
    ero_scatterdata_x=behav_data(s_inds_g(:,group),cond);
    
    if cond==imp.maxconds+1
        ero_scatterdata_y=squeeze(mean(mean(mean(wave_totdata(...
            t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4)) -...
            squeeze(mean(mean(mean(wave_totdata(...
            t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4));
        axis(ero_axes_diff);
    else
        ero_scatterdata_y=squeeze(mean(mean(mean(wave_totdata(...
            t_start:t_end,chan,f_end:f_start,cond,s_inds_g(:,group)),1),3),4)) - ...
            squeeze(mean(mean(mean(wave_totdata(...
            1:scl.t_zero,chan,f_end:f_start,cond,s_inds_g(:,group)),1),3),4));
        axis(ero_axes);
    end

    %p_tfwin(:,chan,freq_range,win)=robustfit(itc_scatterdata_x,itc_scatterdata_y);
    %
    scatter_h(group)=scatter(ero_scatterdata_x,ero_scatterdata_y, scl.g_color{group}); hold on;
    %
    %plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100),'k--'); hold on;
    %plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100)*...
    %    p_tfwin(2,chan,freq_range,win)+p_tfwin(1,chan,freq_range,win), scl.g_color{group}); hold on;
    end
    hold off; grid on;
    end
    end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
tightfig; dragzoom;
end
end
clear_plotassistvars

%% image phase in time-freq at a chosen channel

v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
for group=pp.chosen_g
for chan=pp.chosen_chan(pp.plotn_chan)
figure; subplot_dummy=0;
for cond=1:imp.maxconds
    subplot_dummy=subplot_dummy+1;
    subplot(sp_d(1),sp_d(2),subplot_dummy)
    phase_plot_data=squeeze( angle(mean(wave_evkdata(:,chan,:,cond,s_inds_g(:,group)),5)));
    %phase_plot_data=squeeze( deunwrap(mean(unwrap(angle(wave_evkdata(:,chan,:,cond,s_inds_g(:,group)))),5)));
    contourf(fliplr(contour_wrapedges(phase_plot_data))',pp.n_contour)
    shading flat; colormap(pp.cmap)
    %imagesc(phase_plot_data'); colormap(pp.cmap);
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,cond,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); ylabel('Frequency (Hz)');
    grid on;
    title(['Phase at ',scl.chan_label{chan},' : ',scl.g_label{group},' : ',scl.cond_label{cond}])
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
for splot=1:imp.maxconds;
    subplot(sp_d(1),sp_d(2),splot);
    caxis([-pi pi]); colorbar('YTick',[-pi,0,pi],'YTickLabel',{'0',[char(177),'pi/2'],[char(177),'pi']});
end
tightfig; dragzoom;
end
end
%distFig('s','ext','transpose',true);
clear_plotassistvars

%% plot phase as a topoplot

phase_topo_scale=[-pi pi];
phase_diff_limits=[-pi pi];

for group=pp.chosen_g(pp.plotn_g)
for freq_range=pp.plotn_f
figure; subplot_dummy=0;
%convert scl.freqs to pts
[~,f_indiv]=min(abs(scl.freqs-pp.f_indiv_hz(freq_range)));
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)
        %convert times to points
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        %plot
        subplot_dummy=subplot_dummy+1;
        subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy)
        if cond==imp.maxconds+1
            topo_data=squeeze(angle(mean(mean(mean(wave_evkdata(t_end:t_end,:,f_indiv,pp.cond_diff{1},s_inds_g(:,group)),1),4),5) -...
                mean(mean(mean(wave_evkdata(t_end:t_end,:,f_indiv,pp.cond_diff{2},s_inds_g(:,group)),1),4),5)));
            topoplot(contour_wrapedges(topo_data),chan_locs,'maplimits',[phase_diff_limits(1) phase_diff_limits(2)],'electrodes','off','colormap',pp.cmap);
        else
            topo_data=squeeze(angle(mean(mean(wave_evkdata(t_end:t_end,:,f_indiv,cond,s_inds_g(:,group)),1),5)));
            topoplot(contour_wrapedges(topo_data),chan_locs,'maplimits',[phase_topo_scale(1) phase_topo_scale(2)],'electrodes','off','colormap',pp.cmap);
        end
        title(sprintf('%s, %d - %d ms, %1.1f - %1.1f Hz, %s',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),pp.f_indiv_hz(freq_range),scl.g_label{group}))
    end
    colorbar('YTick',[-pi,0,pi],'YTickLabel',{'0',[char(177),'pi/2'],[char(177),'pi']});
end
tightfig;
end
end
clear_plotassistvars

%% plot phase / ITC at latencies as a wind-rose

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

di_bins=[0.2:0.1:0.7]; %,'di',[0.1:0.1:0.7]); [10:10:50]
di_bins_diff=[-0.3 -0.2 -0.1 0 0.1 0.2 0.3]; %'di',[-0.2:0.05:0.2]); [-20:10:20]

v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_indiv_hz),length(pp.t_start_ms),length(pp.plotn_cond),2);
for group=pp.chosen_g(pp.plotn_g)
for chan=pp.chosen_chan(pp.plotn_chan)
for freq_range=pp.plotn_f;
figure; subplot_dummy=0;
overtitle=sprintf('Phase/ITC Windrose in %1.1f - %1.1f Hz at %s / %s',pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.chan_label{chan},scl.g_label{group});
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
[~,f_indiv]=min(abs(scl.freqs-pp.f_indiv_hz(freq_range)));
for cond=pp.plotn_cond
for win=1:length(pp.t_start_ms)
    [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    sp=subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy);
    if cond==imp.maxconds+1
        phase_data=rad2deg(squeeze(...
            angle(mean(mean(mean(wave_evkdata(t_start:t_end,chan,f_indiv,pp.cond_diff{1},s_inds_g(:,group)),1),3),4)) -...
            angle(mean(mean(mean(wave_evkdata(t_start:t_end,chan,f_indiv,pp.cond_diff{2},s_inds_g(:,group)),1),3),4))));
        intensity_data=squeeze(...
            mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4) - ...
            mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4));
        wind_title=sprintf('%d - %d ms, %1.1f Hz, %s, %s',pp.t_start_ms(win),...
            pp.t_end_ms(win),pp.f_indiv_hz(freq_range),scl.chan_label{chan},...
            scl.g_label{group});
        wind_rose(phase_data,intensity_data,'n',18,'ci',[10 20 30],...
            'lablegend','coherence','cmap',pp.cmap,...
            'quad',3,'parent',sp,'di',di_bins_diff); %'di',[-0.2:0.05:0.2]);
    else
        phase_data=rad2deg(angle(squeeze(mean(mean(wave_evkdata(t_start:t_end, ...
            chan,f_indiv,cond,s_inds_g(:,group)),1),3))));
        intensity_data=squeeze(mean(mean(itcdata(t_start:t_end, ...
            chan,f_end:f_start,cond,s_inds_g(:,group)),1),3));
        wind_title=sprintf('%d - %d ms, %1.1f Hz, %s, %s',pp.t_start_ms(win),...
            pp.t_end_ms(win),pp.f_indiv_hz(freq_range),scl.chan_label{chan},...
            scl.g_label{group});
        wind_rose(phase_data,intensity_data,'n',18,'ci',[10 20 30],...
            'lablegend','coherence','cmap',pp.cmap,...
            'quad',3,'parent',sp,'di',di_bins); %,'di',[0.1:0.1:0.7]);
    end
    v(group,chan,freq_range,win,cond,1)=min(intensity_data);
    v(group,chan,freq_range,win,cond,2)=max(intensity_data);
end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
%tightfig; %dragzoom;
end
end
end
c(1)=min(min(min(min(min(v(:,:,:,:,1:end-1,1))))));
c(2)=max(max(max(max(max(v(:,:,:,:,1:end-1,2))))));
c_diff(1)=min(min(min(min(v(:,:,:,:,end,1)))));
c_diff(2)=max(max(max(max(v(:,:,:,:,end,2)))));
clear_plotassistvars

%% image ITC in time-freq at a chosen channel

pp.figdum=pp.figdum_init;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
for group=pp.chosen_g
for chan=pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(sp_d(1),sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        itc_plot_data=squeeze(mean(mean(itcdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),4),5)-...
            mean(mean(itcdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),4),5));
        %itc_plot_data=squeeze(mean(wavelet_tot(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),5)-...
        %    mean(wavelet_tot(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),5));
    else
        itc_plot_data=squeeze(mean(itcdata(:,chan,:,cond,s_inds_g(:,group)),5));
        %itc_plot_data=squeeze(mean(wavelet_tot(:,chan,:,cond,s_inds_g(:,group)),5));
    end
    contourf(fliplr(itc_plot_data)',pp.n_contour)
    shading flat; colormap(pp.cmap)
    %imagesc(itc_plot_data');
    
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,subplot_dummy,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); ylabel('Frequency (Hz)');
    grid on;
    title(['ITC at ',scl.chan_label{chan},' : ',scl.g_label{group},' : ',scl.cond_label{cond}])
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
end
end
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
for splot=1:length(pp.plotn_cond);
    if splot==length(pp.plotn_cond)
        subplot(sp_d(1),sp_d(2),splot); caxis([c_diff(1) c_diff(2)]);
    else
        subplot(sp_d(1),sp_d(2),splot); caxis([c(1) c(2)]);
    end
    colorbar;
end
tightfig;
end
%distFig('s','ext','transpose',true);
clear_plotassistvars

%% image ITC in time-freq at a chosen channel, with ERP superimposed

pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
for group=pp.chosen_g
for chan=pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(sp_d(1),sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        itc_plot_data=squeeze(mean(mean(itcdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),4),5)-...
            mean(mean(itcdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),4),5));
        %itc_plot_data=squeeze(mean(wavelet_tot(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),5)-...
        %    mean(wavelet_tot(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),5));
        erp_plot_data=mean(mean(erpdata(:,chan,pp.cond_diff{1},s_inds_g(:,group)),3) - ...
            mean(erpdata(:,chan,pp.cond_diff{2},s_inds_g(:,group)),3),4);
    else
        itc_plot_data=squeeze(mean(itcdata(:,chan,:,cond,s_inds_g(:,group)),5));
        erp_plot_data=mean(erpdata(:,chan,cond,s_inds_g(:,group)),4);
        %itc_plot_data=squeeze(mean(wavelet_tot(:,chan,:,cond,s_inds_g(:,group)),5));
    end
    contourf(fliplr(itc_plot_data)',pp.n_contour)
    shading flat; colormap(pp.cmap)
    %imagesc(itc_plot_data');
    hold on;
    %transform ERP to fit?
    plot(erp_plot_data+8,'w');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,subplot_dummy,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); ylabel('Frequency (Hz)');
    grid on;
    title(['ITC at ',scl.chan_label{chan},' : ',scl.g_label{group},' : ',scl.cond_label{cond}])
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--');
end
end
end
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
for splot=1:length(pp.plotn_cond);
    if splot==length(pp.plotn_cond)
        subplot(sp_d(1),sp_d(2),splot); caxis([c_diff(1) c_diff(2)]);
    elseif splot==sp_d(2)
        subplot(sp_d(1),sp_d(2),splot); caxis([c(1) c(2)]);
    end
    colorbar;
end
tightfig;
end
%distFig('s','ext','transpose',true);
clear_plotassistvars

%% image ITC as a scalp plot in a series of time-frequency windows

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

itc_topo_scale=[0.1 0.5];
itc_diff_limits=[-.1 .1];

%
v=zeros(length(pp.plotn_cond),length(pp.chosen_g),length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms),2);
for group=pp.chosen_g(pp.plotn_g)
for freq_range=pp.plotn_f
figure; subplot_dummy=0;
overtitle=sprintf('ITC in %1.1f - %1.1f Hz / %s',pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.g_label{group});
%convert scl.freqs to pts
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)
        %convert times to points
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        %plot
        subplot_dummy=subplot_dummy+1;
        subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy)
        if cond==imp.maxconds+1
            topo_data=squeeze(mean(mean(mean(mean(itcdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4),5)) - ...
                squeeze(mean(mean(mean(mean(itcdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4),5));
            topoplot(topo_data,chan_locs,'maplimits',[itc_diff_limits(1) itc_diff_limits(2)],'electrodes','off','colormap',pp.cmap);
        else
            topo_data=squeeze(mean(mean(mean(itcdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5));
            topoplot(topo_data,chan_locs,'maplimits',[itc_topo_scale(1) itc_topo_scale(2)],'electrodes','off','colormap',pp.cmap);
        end
        v(cond,group,freq_range,win,1) = min(topo_data); v(cond,group,freq_range,win,2) = max(topo_data);
        %title(sprintf('%s, %d - %d ms, %1.1f - %1.1f Hz, %s',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.g_label{group}))
    end
    colorbar;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
%tightfig;
%dragzoom;
end
end
c(1)=min(min(min(min(v(1:end-1,:,:,:,1))))); c(2)=max(max(max(max(v(1:end-1,:,:,:,2)))));
c_diff(1)=min(min(min(v(end,:,:,:,1)))); c_diff(2)=max(max(max(v(end,:,:,:,2))));
%c(2)=itc_topo_scale;%c(2) override
%for splot=1:subplot_dummy; subplot(imp.maxconds,length(pp.t_start_ms),splot); caxis([c(1) c(2)]); end
clear_plotassistvars

%% ITC - create condition bar plots with error bars

sp_rowlabel=make_freqlabels(pp.f_start_hz(pp.plotn_f),pp.f_end_hz(pp.plotn_f));
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms)];

width=[];
bw_xlabel=[];
bw_ylabel='ITC';
bw_colormap=bone; %pmkmp(length(barvalues),'cubicl');
gridstatus='y';
error_sides=1;
legend_type='plot'; %'plot'
%legend_type=[];
bar_glabel={scl.g_label{pp.chosen_g(pp.plotn_g)}};
%bar_condlabel={'Go','NoGo'};
bar_condlabel={scl.cond_label{pp.plotn_cond(1:end-1)}};
%bar_condlabel=[];

%figures are chans, columns are time windows, rows are frequency bands
for chan=pp.chosen_chan(pp.plotn_chan)
subplot_dummy=0;
figure;
overtitle=sprintf('BarplotWEB of ITC for %s',scl.chan_label{chan});
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    for win=1:length(pp.t_start_ms)
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        subplot_dummy=subplot_dummy+1;
        subplot(length(pp.plotn_f),pp.maxwin,subplot_dummy); %length(pp.f_start_hz(pp.plotn_f)),pp.maxwin,subplot_dummy);
        %
        gdum=0;
        for group=pp.chosen_g(pp.plotn_g)
        gdum=gdum+1;
        ranova_data{gdum}=squeeze(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.plotn_cond(1:end-1),s_inds_g(:,group)),1),3))';
        bar_data(group,:)=squeeze(mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.plotn_cond(1:end-1),s_inds_g(:,group)),1),3),5));
        bar_data_se(group,:)=std(squeeze(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.plotn_cond(1:end-1),s_inds_g(:,group)),1),3)),0,2)/sqrt(sum(s_inds_g(:,group)));
        %
        end
        %bar_title=sprintf('%s, %1.1f - %1.1f, %d - %d ms',scl.chan_label{chan},...
        %    pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),pp.t_start_ms(win),pp.t_end_ms(win));
        bar_title=[];
        bar_h=barweb(bar_data(pp.chosen_g,:),bar_data_se(pp.chosen_g,:),width,bar_glabel,...
            bar_title,bw_xlabel,bw_ylabel,bw_colormap,gridstatus,bar_condlabel,error_sides,legend_type);
        axis 'auto y'
        axis([0 4 0.2 0.6])
        [p,ranova_table]=anova_rm(ranova_data,'off');
        text(2.5,0.32,sprintf('main p=%1.3f',p(1)))
        text(2.5,0.25,sprintf('group p=%1.3f',p(2)))
        %text(2.5,0.18,sprintf('int p=%1.3f',p(3)))
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
end
%linkaxes
tightfig;
clear_plotassistvars

%% scatter ITC with RT measures

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='RT Variance (ms)'; %Median, variance, std also interesting
y_plotlabel='ITC';

behav_data=stdRT;

itc_axes=[min(min(behav_data)) max(max(behav_data)) 0.1 0.8]; %max(max(behav_data))
itc_axes_diff=[min(min(behav_data)) max(max(behav_data)) -2 5];

p_tfwin=zeros(2,length(pp.plotn_cond),length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_start_hz(pp.plotn_f)),pp.maxwin);
for chan=pp.chosen_chan(pp.plotn_chan)
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    figure; subplot_dummy=0;
    overtitle=sprintf('%s in %1.1f - %1.1f Hz', ...
        scl.chan_label{chan},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
    for cond=pp.plotn_cond
    for win=1:pp.maxwin
    [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
    %    
    for group=pp.chosen_g(pp.plotn_g)
        
    itc_scatterdata_x=behav_data(s_inds_g(:,group),cond);
    
    if cond==imp.maxconds+1
        itc_scatterdata_y=squeeze(mean(mean(mean(itcdata(...
            t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4)) -...
            squeeze(mean(mean(mean(itcdata(...
            t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4));
        axis(itc_axes_diff);
    else
        itc_scatterdata_y=squeeze(mean(mean(mean(itcdata(...
            t_start:t_end,chan,f_end:f_start,cond,s_inds_g(:,group)),1),3),4));
        axis(itc_axes);
    end

    %p_tfwin(:,cond,chan,freq_range,win)=robustfit(itc_scatterdata_x,itc_scatterdata_y);
    p_tfwin(:,cond,chan,freq_range,win)=regress(itc_scatterdata_y,[ones(length(itc_scatterdata_x),1) itc_scatterdata_x]);
    %
    scatter_h(group)=scatter(itc_scatterdata_x,itc_scatterdata_y, scl.g_color{group}); hold on;
    %
    %plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100),'k--'); hold on;
    plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100)*...
        p_tfwin(2,cond,chan,freq_range,win)+p_tfwin(1,cond,chan,freq_range,win), scl.g_color{group}); hold on;
    end
    hold off; grid on;
    end
    end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
tightfig; dragzoom;
end
end
clear_plotassistvars

%% SCATTER ITC

sp_rowlabel=make_freqlabels(pp.f_start_hz(pp.plotn_f),pp.f_end_hz(pp.plotn_f));
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=[scl.cond_label{pp.cond_diff{1}}];
y_plotlabel=[scl.cond_label{pp.cond_diff{2}}];
subplot_dims=[length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms)];

itc_axes=[0 1 0 1];

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_start_hz(pp.plotn_f)),pp.maxwin);
for chan=pp.chosen_chan(pp.plotn_chan)
    figure; subplot_dummy=0;
    overtitle=sprintf('ITC for %s',scl.chan_label{chan});
    for freq_range=pp.plotn_f
        [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
        [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
        for win=1:length(pp.t_start_ms)
            subplot_dummy=subplot_dummy+1;
            subplot(length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms),subplot_dummy);
            %
            [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
            [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
            %
            for group=pp.chosen_g(pp.plotn_g)
                itc_scatterdata_x=squeeze(mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4));
                itc_scatterdata_y=squeeze(mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4));
                %itc_scatterdata_x=squeeze(mean(mean(itcdata(t_start:t_end,chan,18:20,pp.cond_diff{1},s_inds_g(:,group)),1),3)-mean(mean(itcdata(t_start:t_end,chan,18:20,pp.cond_diff{2},s_inds_g(:,group)),1),3));
                %itc_scatterdata_y=squeeze(mean(mean(itcdata(t_start:t_end,chan,12:14,pp.cond_diff{1},s_inds_g(:,group)),1),3)-mean(mean(itcdata(t_start:t_end,chan,12:14,pp.cond_diff{2},s_inds_g(:,group)),1),3));
                %p_tfwin(:,chan,freq_range,win)=polyfit(itc_scatterdata_x,itc_scatterdata_y,1);
                p_tfwin(:,chan,freq_range,win)=robustfit(itc_scatterdata_x,itc_scatterdata_y);
                %
                scatter_h(group)=scatter(itc_scatterdata_x,itc_scatterdata_y, scl.g_color{group}); hold on;
                %scatter(s_demogs(s_inds_g(:,group),:).EEG_Age,itc_scatterdata_x, scl.g_color{group-4}); hold on;
                %scatter(s_demogs(s_inds_g(:,group),:).EEG_Age,itc_scatterdata_y, scl.g_color{group-2}); hold on;
                %scatter(itc_scatterdata_x+itc_scatterdata_y,itc_scatterdata_x-itc_scatterdata_y); hold on; axis(itc_axes); grid on;
                plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100),'k--'); hold on;
                plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100)*p_tfwin(2,chan,freq_range,win)+p_tfwin(1,chan,freq_range,win), scl.g_color{group}); hold on;
            end
            hold off; axis(itc_axes); grid on;
            %title(sprintf('%d - %d ms, %1.1f - %1.1f Hz, %s',pp.t_start_ms(win),pp.t_end_ms(win),pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.chan_label{chan}));
            %xlabel([scl.cond_label{pp.cond_diff{1}}]); ylabel([scl.cond_label{pp.cond_diff{2}}])
            %xlabel([scl.cond_label{pp.cond_diff{1}},'+',scl.cond_label{pp.cond_diff{2}}]); ylabel([scl.cond_label{pp.cond_diff{1}},'-',scl.cond_label{pp.cond_diff{2}}]);
            %xlabel('Lo Theta Diff'); ylabel('Hi Theta Diff')
            %legend(scatter_h(pp.chosen_g),scl.g_label{pp.chosen_g})
        end
    end
    %tightfig; dragzoom;
    adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
end
clear_plotassistvars

%% SCATTER ERO WITH AGE

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='Age (yrs)';
y_plotlabel='ERO Power';

ero_axes=[min(s_demogs.age_eeg) max(s_demogs.age_eeg) 40 90];
ero_axes_diff=[min(s_demogs.age_eeg) max(s_demogs.age_eeg) -10 10];

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_start_hz(pp.plotn_f)),pp.maxwin);
for chan=pp.chosen_chan(pp.plotn_chan)
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    figure; subplot_dummy=0;
    overtitle=sprintf('ERO Power over age at %s in %1.1f - %1.1f Hz', ...
        scl.chan_label{chan},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
    for cond=pp.plotn_cond
    for win=1:pp.maxwin
    [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
    %    
    for group=pp.chosen_g(pp.plotn_g)
    ero_scatterdata_x=s_demogs.age_eeg(s_inds_g(:,group));
    if cond==imp.maxconds+1
        ero_scatterdata_y=squeeze(mean(mean(mean(wave_totdata(...
            t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4)) -...
            squeeze(mean(mean(mean(wave_totdata(...
            t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4));
        axis(ero_axes_diff);
    else
        ero_scatterdata_y=squeeze(mean(mean(mean(wave_totdata(...
            t_start:t_end,chan,f_end:f_start,cond,s_inds_g(:,group)),1),3),4));
        axis(ero_axes);
    end

    %p_tfwin(:,chan,freq_range,win)=robustfit(itc_scatterdata_x,itc_scatterdata_y);
    %
    scatter_h(group)=scatter(ero_scatterdata_x,ero_scatterdata_y, scl.g_color{group}); hold on;
    %
    %plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100),'k--'); hold on;
    %plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100)*...
    %    p_tfwin(2,chan,freq_range,win)+p_tfwin(1,chan,freq_range,win), scl.g_color{group}); hold on;
    end
    hold off; grid on;
    end
    end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
%tightfig; dragzoom;
end
end
clear_plotassistvars


%% SCATTER ITC WITH AGE

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='Age (yrs)';
y_plotlabel='ITC';

itc_axes=[min(s_demogs.age_eeg) max(s_demogs.age_eeg) 0 0.8];
itc_axes_diff=[min(s_demogs.age_eeg) max(s_demogs.age_eeg) -0.5 0.5];

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_start_hz(pp.plotn_f)),pp.maxwin);
for chan=pp.chosen_chan(pp.plotn_chan)
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    figure; subplot_dummy=0;
    overtitle=sprintf('ITC over age at %s in %1.1f - %1.1f Hz', ...
        scl.chan_label{chan},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
    for cond=pp.plotn_cond
    for win=1:pp.maxwin
    [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
    %    
    for group=pp.chosen_g(pp.plotn_g)
    itc_scatterdata_x=s_demogs.age_eeg(s_inds_g(:,group));
    if cond==imp.maxconds+1
        itc_scatterdata_y=squeeze(mean(mean(mean(itcdata(...
            t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4)) -...
            squeeze(mean(mean(mean(itcdata(...
            t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4));
        axis(itc_axes_diff);
    else
        itc_scatterdata_y=squeeze(mean(mean(mean(itcdata(...
            t_start:t_end,chan,f_end:f_start,cond,s_inds_g(:,group)),1),3),4));
        axis(itc_axes);
    end

    %p_tfwin(:,chan,freq_range,win)=robustfit(itc_scatterdata_x,itc_scatterdata_y);
    %
    scatter_h(group)=scatter(itc_scatterdata_x,itc_scatterdata_y, scl.g_color{group}); hold on;
    %
    %plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100),'k--'); hold on;
    %plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100)*...
    %    p_tfwin(2,chan,freq_range,win)+p_tfwin(1,chan,freq_range,win), scl.g_color{group}); hold on;
    end
    hold off;  grid on;
    end
    end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
tightfig; dragzoom;
end
end
clear_plotassistvars

%% HISTOGRAM ITC DIFFERENCE
% PAIRS

sp_rowlabel=make_freqlabels(pp.f_start_hz(pp.plotn_f),pp.f_end_hz(pp.plotn_f));
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.t_start_ms),pp.maxwin);
for group=pp.chosen_g(pp.plotn_g)
for chan=pp.chosen_chan(pp.plotn_chan)
    figure; subplot_dummy=0;
    overtitle=sprintf('ITC for %s at %s in %s',scl.cond_label{end},scl.chan_label{chan},scl.g_label{group});
    for freq_range=pp.plotn_f
        [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
        [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
        for win=1:pp.maxwin
            subplot_dummy=subplot_dummy+1;
            [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
            [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
            %
            itc_histdata=squeeze(mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4) - ...
                mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4));
            itc_histdata_se=std(squeeze(mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4) - ...
                mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4))); %/sqrt(sum(s_inds_g(:,group)));
            %
            subplot(length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms),subplot_dummy);
            hist(itc_histdata,pp.hist_nbin); axis(pp.hist_ax); hold on;
            h=findobj(gca,'Type','patch'); set(h,'FaceColor',[.7 .7 .7],'EdgeColor','w');
            plot(ones(1,100)*mean(itc_histdata),linspace(pp.hist_ax(3),pp.hist_ax(4),100),'k--'); hold on;
            plot(linspace(mean(itc_histdata)-itc_histdata_se,mean(itc_histdata)+itc_histdata_se,100), ...
                ones(1,100)*pp.hist_ax(4)/20,'k','LineWidth',2); hold off;
            %[f,x]=ksdensity(itc_histdata); plot(x,f);
            [~,p,~,stat]=ttest(zeros(sum(s_inds_g(:,group)),1),itc_histdata);
            z=mean(itc_histdata)/itc_histdata_se;
            text(pp.p_loc(1),pp.p_loc(2),sprintf('z=%1.1f, t=%1.1f, -log(p)=%1.1f',z,stat.tstat,-log10(p)));
            %title(sprintf('%d - %d ms, %1.1f - %1.1f Hz, %s, %s',pp.t_start_ms(win),pp.t_end_ms(win),pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.chan_label{chan},scl.g_label{group}));
            %xlabel([scl.cond_label{pp.cond_diff{1}},'-',scl.cond_label{pp.cond_diff{2}}])
        end
    end
    adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
    tightfig; dragzoom;
end
end
clear_plotassistvars

%% STACKED HISTOGRAM ITC DIFFERENCE

sp_rowlabel=make_freqlabels(pp.f_start_hz(pp.plotn_f),pp.f_end_hz(pp.plotn_f));
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.t_start_ms),pp.maxwin);

for chan=pp.chosen_chan(pp.plotn_chan)
figure; subplot_dummy=0;
overtitle=sprintf('ITC for %s at %s',scl.cond_label{end},scl.chan_label{chan});
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    for win=1:pp.maxwin
        subplot_dummy=subplot_dummy+1;
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        %
        for group=pp.chosen_g(pp.plotn_g)
            itc_histdata=squeeze(mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4) - ...
                mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4));
            itc_histdata_se=std(squeeze(mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4) - ...
                mean(mean(mean(itcdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4))); %/sqrt(sum(s_inds_g(:,group)));
            %
            [nels(:,group),cents(:,group)]=hist(itc_histdata,pp.hist_nbin);
        end
        subplot(length(pp.f_start_hz(pp.plotn_f)),pp.maxwin,subplot_dummy);
        bar(cents(:,pp.chosen_g),nels(:,pp.chosen_g),'stacked');
        axis(pp.hist_ax); hold on;
        h=findobj(gca,'Type','patch');
        set(h(1),'FaceColor',[.7 .7 .7],'EdgeColor','w');
        set(h(2),'FaceColor',[.5 .5 .5],'EdgeColor','w');
        plot(ones(1,100)*mean(itc_histdata),linspace(pp.hist_ax(3),pp.hist_ax(4),100),'k--'); hold on;
        plot(linspace(mean(itc_histdata)-itc_histdata_se,mean(itc_histdata)+itc_histdata_se,100), ...
            ones(1,100)*pp.hist_ax(4)/20,'k','LineWidth',2); hold off;
        %[f,x]=ksdensity(itc_histdata); plot(x,f);
        [~,p,~,stat]=ttest(zeros(sum(s_inds_g(:,group)),1),itc_histdata);
        z=mean(itc_histdata)/itc_histdata_se;
        text(pp.p_loc(1),pp.p_loc(2),sprintf('z=%1.1f, t=%1.1f, -log(p)=%1.1f',z,stat.tstat,-log10(p)));
        %title(sprintf('%d - %d ms, %1.1f - %1.1f Hz, %s',pp.t_start_ms(win),pp.t_end_ms(win),pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.chan_label{chan}));
        xlabel([scl.cond_label{pp.cond_diff{1}},'-',scl.cond_label{pp.cond_diff{2}}])
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
tightfig; dragzoom;
end
clear_plotassistvars

%% KSDENSITY ITC

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='ITC';
y_plotlabel='proportion of subjects';

for chan=pp.chosen_chan(pp.plotn_chan)
for freq_range=pp.plotn_f
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
figure; subplot_dummy=0;
overtitle=sprintf('Kernel Density of ITC at %s for %1.1f - %1.1f Hz', ...
    scl.chan_label{chan},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
for cond=pp.plotn_cond
for win=1:pp.maxwin
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
    [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
    %
    for group=pp.chosen_g(pp.plotn_g)
        if cond==imp.maxconds+1
            itc_histdata=squeeze(mean(mean(mean(itcdata(t_start:t_end,chan, ...
            f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4)) - ...
            squeeze(mean(mean(mean(itcdata(t_start:t_end,chan, ...
            f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4));
            itc_histdata_se=std(squeeze(mean(mean(mean(itcdata(t_start:t_end, ...
            chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4)) -...
            squeeze(mean(mean(mean(itcdata(t_start:t_end, ...
            chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4)));
        else
            itc_histdata=squeeze(mean(mean(mean(itcdata(t_start:t_end,chan, ...
            f_end:f_start,cond,s_inds_g(:,group)),1),3),4));
            itc_histdata_se=std(squeeze(mean(mean(mean(itcdata(t_start:t_end, ...
            chan,f_end:f_start,cond,s_inds_g(:,group)),1),3),4)));
        end
        %
        [f,xi]=ksdensity(itc_histdata);
        f=f/max(abs(f))*pp.hist_ax(4);
        %
        plot(xi,f,scl.g_color{group}); hold on;
        plot(ones(1,100)*mean(itc_histdata),linspace(pp.hist_ax(3),pp.hist_ax(4),100),[scl.g_color{group},'--']);
        hold on; plot(linspace(mean(itc_histdata)-itc_histdata_se,mean(itc_histdata)+itc_histdata_se,100), ...
        ones(1,100)*pp.hist_ax(4)/20+(group-10),scl.g_color{group},'LineWidth',2); hold on;
    end
    axis(pp.hist_ax); hold on;
    %[~,p,~,stat]=ttest(itc_histdata,itc_histdata);
    %z=mean(itc_histdata)/itc_histdata_se;
    %text(pp.p_loc(1),pp.p_loc(2),sprintf('z=%1.1f, t=%1.1f, -log(p)=%1.1f',z,stat.tstat,-log10(p)));
end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
tightfig; dragzoom;
end
end
clear_plotassistvars

%% CORRELATION STRUCTURE OF FREQUENCY BANDS ACROSS TIME, ERPCOH



%% CORRELATION STRUCTURE OF FREQUENCY BANDS ACROSS TIME FOR ITC

for group=pp.chosen_g
v=zeros(length(pp.plotn_cond),pp.maxwin,2);
for chan=pp.chosen_chan(pp.plotn_chan)
figure; subplot_dummy=0;
for cond=pp.plotn_cond
    for win=1:pp.maxwin
        subplot_dummy=subplot_dummy+1;
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        %
        if cond==imp.maxconds+1
            itc_freqcorr_data=squeeze(mean(mean(itcdata(t_start:t_end,chan,:,pp.cond_diff{1},s_inds_g(:,group)),1),4)-...
                mean(mean(itcdata(t_start:t_end,chan,:,pp.cond_diff{2},s_inds_g(:,group)),1),4))';
            %itc_freqcorr_data=squeeze(mean(mean(itcdata(t_start:t_end,:,:,pp.cond_diff{1},s_inds_g(:,group)),1),2)-...
             %    mean(mean(itcdata(t_start:t_end,:,:,pp.cond_diff{2},s_inds_g(:,group)),1),2))';
        else
            itc_freqcorr_data=squeeze(mean(itcdata(t_start:t_end,chan,:,cond,s_inds_g(:,group)),1))';
            %itc_freqcorr_data=squeeze(mean(mean(itcdata(t_start:t_end,:,:,cond,s_inds_g(:,group)),1),2))';
        end
        itc_freqcorr=corr(itc_freqcorr_data);
        itc_freqcorr=dac(itc_freqcorr);
        %itc_freqcorr(itc_freqcorr < 0) = 0;
        %itc_freqcorr(itc_freqcorr > 0) = 0;
        itc_freqcorr=tril(itc_freqcorr,-1);
        %
        subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
        %imagesc(itc_freqcorr); %caxis(rho_limits);
        %axis xy;
        contourf(itc_freqcorr); %caxis(rho_limits);
        colormap(pp.cmap)
        v(cond,win,:)=caxis;
        title(sprintf('%s, %d - %d ms, %s',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),scl.chan_label{chan}));
        %title(sprintf('%s, %d - %d ms',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win)));
        xlabel('Freq. (Hz)'); ylabel('Freq. (Hz)');
    end
end
c(1)=min(min(v(:,:,1))); c(2)=max(max(v(:,:,2)));
for splot=1:subplot_dummy
    subplot(length(pp.plotn_cond),pp.maxwin,splot); caxis([c(1) c(2)]);
end
tightfig; %dragzoom;
end
end
%distFig;
clear_plotassistvars

%% plot "inter-subject phase consistency" as lines

sp_rowlabel=make_freqlabels(pp.f_start_hz(pp.plotn_f),pp.f_end_hz(pp.plotn_f));
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='Time Window';
y_plotlabel='ISPC';

for group=pp.chosen_g(pp.plotn_g)
for chan=pp.chosen_chan(pp.plotn_chan)
figure; subplot_dummy=0;
overtitle=sprintf('ISPC at %s',scl.chan_label{chan});
for freq_range=pp.plotn_f
    [~,f_indiv]=min(abs(scl.freqs-pp.f_indiv_hz(freq_range)));
    for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.f_indiv_hz(pp.plotn_f)),length(pp.plotn_cond),subplot_dummy);
    if cond==imp.maxconds+1
        ispc_data=mean(abs(sum(wave_evkdata(:,chan,f_indiv,pp.cond_diff{1},s_inds_g(:,group)),5)),4)./ ...
            mean(sum(abs(wave_evkdata(:,chan,f_indiv,pp.cond_diff{1},s_inds_g(:,group))),5),4) - ...
            mean(abs(sum(wave_evkdata(:,chan,f_indiv,pp.cond_diff{2},s_inds_g(:,group)),5)),4)./ ...
            mean(sum(abs(wave_evkdata(:,chan,f_indiv,pp.cond_diff{2},s_inds_g(:,group))),5),4);
        plot(ispc_data,'r'); hold on;
        ispc_data=mean(abs(sum(exp(1i*angle(wave_evkdata(:,chan,f_indiv,pp.cond_diff{1},s_inds_g(:,group)))),5)),4)./ ...
            mean(sum(abs(exp(1i*angle(wave_evkdata(:,chan,f_indiv,pp.cond_diff{1},s_inds_g(:,group))))),5),4) - ...
            mean(abs(sum(exp(1i*angle(wave_evkdata(:,chan,f_indiv,pp.cond_diff{2},s_inds_g(:,group)))),5)),4)./ ...
            mean(sum(abs(exp(1i*angle(wave_evkdata(:,chan,f_indiv,pp.cond_diff{2},s_inds_g(:,group))))),5),4);
        plot(ispc_data,'g'); hold on;
    else
        ispc_data=abs(sum(wave_evkdata(:,chan,f_indiv,cond,s_inds_g(:,group)),5))./ ...
            sum(abs(wave_evkdata(:,chan,f_indiv,cond,s_inds_g(:,group))),5);
        plot(ispc_data,'r'); hold on;
        ispc_data=abs(sum(exp(1i*angle(wave_evkdata(:,chan,f_indiv,cond,s_inds_g(:,group)))),5))./ ...
            sum(abs(exp(1i*angle(wave_evkdata(:,chan,f_indiv,cond,s_inds_g(:,group))))),5);
        plot(ispc_data,'g'); hold on;
    end
    if cond==imp.maxconds+1
    grid on; axis([scl.t_start scl.t_end pp.chosen_coh_limit_diff(1) pp.chosen_coh_limit_diff(2)]);
    else
    grid on; axis([scl.t_start scl.t_end pp.chosen_coh_limit(1) pp.chosen_coh_limit(2)]);
    end
    plot(ones(10,1)*scl.t_zero,linspace(-1,1,10),'k--'); hold off;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
    %title(sprintf('%s, %s, %1.1f Hz,%s',scl.chan_label{chan},scl.cond_label{cond},pp.f_indiv_hz(freq_range),scl.g_label{group}))
    %clickableLegend('weighted','non-weighted')
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
end
end
clear_plotassistvars

%% image "inter-subject phase consistency" in time-freq at a chosen channel

pp.figdum=pp.figdum_init;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
for group=pp.chosen_g
for chan=pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(sp_d(1)+1,sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        ispc_plot_data=squeeze( mean(abs(sum(wave_evkdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),5))./ ...
            sum(abs(wave_evkdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group))),5),4) ) - ...
            squeeze( mean(abs(sum(wave_evkdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),5))./ ...
            sum(abs(wave_evkdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group))),5),4) );
    else
        ispc_plot_data=squeeze( abs(sum(wave_evkdata(:,chan,:,cond,s_inds_g(:,group)),5))./ ...
            sum(abs(wave_evkdata(:,chan,:,cond,s_inds_g(:,group))),5) );
    end
    contourf(fliplr(ispc_plot_data)',pp.n_contour)
    shading flat; colormap(pp.cmap)
    %imagesc(ispc_plot_data');
    
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,cond,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); ylabel('Frequency (Hz)');
    grid on;
    title(['ISPC at ',scl.chan_label{chan},' : ',scl.g_label{group},' : ',scl.cond_label{cond}])
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
end
end
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
for splot=pp.plotn_cond;
    if splot==imp.maxconds+1
        subplot(sp_d(1)+1,sp_d(2),splot); caxis([c_diff(1) c_diff(2)]);
    else
        subplot(sp_d(1)+1,sp_d(2),splot); caxis([c(1) c(2)]);
    end
    colorbar;
end
tightfig;
end
%distFig('s','ext','transpose',true);
clear_plotassistvars

%% plot "inter-subject phase consistency" as a topoplot over time

ispc_topo_scale=[0 0.9];
ispc_diff_limits=[-.5 .5];

v=zeros(length(pp.plotn_cond),length(pp.f_indiv_hz),length(pp.t_start_ms),length(pp.chosen_g),2);
for group=pp.chosen_g(pp.plotn_g)
for freq_range=pp.plotn_f
    figure; subplot_dummy=0;
    [~,f_indiv]=min(abs(scl.freqs-pp.f_indiv_hz(freq_range)));
    for cond=pp.plotn_cond
    for win=1:pp.maxwin
        subplot_dummy=subplot_dummy+1;
        subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
        %
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        %
        if cond==imp.maxconds+1
            ispc_data=squeeze(mean(mean(mean(abs(sum(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,pp.cond_diff{1},s_inds_g(:,group)),5))./ ...
                sum(abs(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,pp.cond_diff{1},s_inds_g(:,group))),5),4),3),1)) - ...
                squeeze(mean(mean(mean(abs(sum(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,pp.cond_diff{2},s_inds_g(:,group)),5))./ ...
                sum(abs(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,pp.cond_diff{2},s_inds_g(:,group))),5),4),3),1));
            %ispc_data=squeeze(mean(mean(mean(abs(sum(exp(1i*angle(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,pp.cond_diff{1},s_inds_g(:,9)))),5))./ ...
            %    sum(abs(exp(1i*angle(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,pp.cond_diff{1},s_inds_g(:,9))))),5),4),3),1)) - ...
            %    squeeze(mean(mean(mean(abs(sum(exp(1i*angle(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,pp.cond_diff{2},s_inds_g(:,9)))),5))./ ...
            %    sum(abs(exp(1i*angle(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,pp.cond_diff{2},s_inds_g(:,9))))),5),4),3),1));
            topoplot(ispc_data,chan_locs,'maplimits',[ispc_diff_limits(1) ispc_diff_limits(2)],'electrodes','off','colormap',pp.cmap);
        else
            ispc_data=squeeze(mean(mean(abs(sum(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,cond,s_inds_g(:,group)),5))./ ...
                sum(abs(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,cond,s_inds_g(:,group))),5),3),1));
            %ispc_data=squeeze(mean(mean(abs(sum(exp(1i*angle(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,cond,s_inds_g(:,9)))),5))./ ...
            %    sum(abs(exp(1i*angle(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,cond,s_inds_g(:,9))))),5),3),1));
            topoplot(ispc_data,chan_locs,'maplimits',[ispc_topo_scale(1) ispc_topo_scale(2)],'electrodes','off','colormap',pp.cmap);
        end
        v(cond,freq_range,win,group,1) = min(ispc_data); v(cond,freq_range,win,group,2) = max(ispc_data);
        title(sprintf('%s, %d - %d ms, %1.1f Hz, %s',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),pp.f_indiv_hz(freq_range),scl.g_label{group}));
    end
    colorbar;
    end
end
end
c(1)=min(min(min(min(v(1:end-1,:,:,pp.chosen_g,1))))); c(2)=max(max(max(max(v(1:end-1,:,:,pp.chosen_g,2)))));
c_diff(1)=min(min(min(v(end,:,:,pp.chosen_g,1)))); c_diff(2)=max(max(max(v(end,:,:,pp.chosen_g,2))));
clear_plotassistvars

%% SCATTER DIFFERENCE FROM MEAN PHASE VS. ITC

phase_axes=[-pi pi 0 1];

sp_rowlabel={scl.cond_label{1:imp.maxconds}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='difference from mean phase';
y_plotlabel='ITC';

p_phase=zeros(2,length(pp.plotn_cond),pp.maxwin);
for chan=pp.chosen_chan(pp.plotn_chan)
for freq=pp.plotn_f
[~,f_start]=min(abs(scl.freqs-pp.f_indiv_hz(freq)));
[~,f_end]=min(abs(scl.freqs-pp.f_indiv_hz(freq)));
[~,f_indiv]=min(abs(scl.freqs-pp.f_indiv_hz(freq)));
figure; subplot_dummy=0;
overtitle=sprintf('Difference from Mean Phase vs. ITC at %s for %1.1f Hz',scl.chan_label{chan},pp.f_indiv_hz(freq));
for cond=1:imp.maxconds
    for win=1:pp.maxwin
        subplot_dummy=subplot_dummy+1;
        subplot(imp.maxconds,pp.maxwin,subplot_dummy);
        %
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        %
        for group=pp.chosen_g(pp.plotn_g)
        if cond==imp.maxconds+1
            scatterdata_x=(squeeze(mean(angle(wave_evkdata(t_start,chan,f_indiv,pp.cond_diff{1},s_inds_g(:,group))),4)) - ...
                mean(mean(angle(wave_evkdata(t_start,chan,f_indiv,pp.cond_diff{1},s_inds_g(:,scl.g_all))),4),5)) - ...
                (squeeze(mean(angle(wave_evkdata(t_start,chan,f_indiv,pp.cond_diff{2},s_inds_g(:,group))),4)) - ...
                mean(mean(angle(wave_evkdata(t_start,chan,f_indiv,pp.cond_diff{2},s_inds_g(:,scl.g_all))),4),5));
            scatterdata_y=squeeze(mean(itcdata(t_start,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),4) - ...
                mean(itcdata(t_start,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),4));
        else
            scatterdata_x=deunwrap(squeeze(unwrap(angle(wave_evkdata(t_start,chan,f_indiv,cond,s_inds_g(:,group))))) - ...
                angle(mean(wave_evkdata(t_start,chan,f_indiv,cond,s_inds_g(:,scl.g_all)),5)));
            %scatterdata_x=squeeze(angle(wave_evkdata(t_start,chan,f_indiv,cond,s_inds_g(:,group)))) - ...
            %    mean(angle(wave_evkdata(t_start,chan,f_indiv,cond,s_inds_g(:,group))),5);
            scatterdata_y=squeeze(itcdata(t_start,chan,f_end:f_start,cond,s_inds_g(:,group)));
        end
        %
        p_phase(:,cond,win)=polyfit(scatterdata_x,scatterdata_y,1);
        %
        scatter(scatterdata_x,scatterdata_y,scl.g_color{group}); hold on;
        %
        %plot(linspace(phase_axes(1),phase_axes(2),100),linspace(phase_axes(3),phase_axes(4),100),'k--'); hold on;
        %plot(linspace(phase_axes(1),phase_axes(2),100), ...
        %    linspace(phase_axes(3),phase_axes(4),100)*p_phase(1,cond,win)+p_phase(2,cond,win),scl.g_color{group});
        end
        hold off; axis(phase_axes); grid on;        
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
tightfig;
end
end
clear_plotassistvars

%% SCATTER DIFFERENCE FROM MEAN PHASE at t1 VS. DIFFERENCE FROM MEAN PHASE at t2

phase_axes=[-pi pi -pi pi];

p_phase=zeros(2,length(pp.plotn_cond),pp.maxwin);
for chan=pp.chosen_chan(pp.plotn_chan)
for freq=1:length(pp.f_indiv_hz)
[~,f_indiv]=min(abs(scl.freqs-pp.f_indiv_hz(freq)));
figure; subplot_dummy=0;
for cond=1:imp.maxconds
    for win=1:pp.maxwin
        subplot_dummy=subplot_dummy+1;
        subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
        %
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        %
        for group=pp.chosen_g(pp.plotn_g)
        scatterdata_x=deunwrap(squeeze(angle(wave_evkdata(t_end,chan,f_indiv,cond,s_inds_g(:,group)))) - ...
            angle(mean(wave_evkdata(t_end,chan,f_indiv,cond,s_inds_g(:,scl.g_all)),5)));
        scatterdata_y=deunwrap(squeeze(angle(wave_evkdata(t_start,chan,f_indiv,cond,s_inds_g(:,group)))) - ...
            angle(mean(wave_evkdata(t_start,chan,f_indiv,cond,s_inds_g(:,scl.g_all)),5)));
        %p_phase(:,cond,win)=polyfit(scatterdata_x,scatterdata_y,1);
        p_phase(:,cond,win)=robustfit(scatterdata_x,scatterdata_y);
        %
        scatter(scatterdata_x,scatterdata_y,scl.g_color{group}); hold on;
        %
        plot(linspace(phase_axes(1),phase_axes(2),100),linspace(phase_axes(3),phase_axes(4),100),'k--'); hold on;
        plot(linspace(phase_axes(1),phase_axes(2),100), ...
            linspace(phase_axes(3),phase_axes(4),100)*p_phase(2,cond,win)+p_phase(1,cond,win),scl.g_color{group});
        end
        hold off; axis(phase_axes); grid on;
        title(sprintf('%1.1f Hz, %s',pp.f_indiv_hz(freq),scl.chan_label{chan}));
        xlabel(sprintf('diff at t=%d',pp.t_end_ms(win)));
        ylabel(sprintf('diff at t=%d',pp.t_start_ms(win)));
    end
end
tightfig; dragzoom;
end
end
clear_plotassistvars