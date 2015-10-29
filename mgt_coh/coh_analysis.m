% SCATTER ERPCOH
% PAIRS

scatter_limits=[0 1];

sp_rowlabel={scl.f_nlabel{pp.chosen_freq}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=scl.cond_label{pp.cond_diff{1}};
y_plotlabel=scl.cond_label{pp.cond_diff{2}};

%s_high_coh=zeros(length(pp.chosen_p),length(pp.chosen_freq),pp.maxwin,10);
p_tfwin=zeros(2,length(pp.chosen_p),length(pp.chosen_freq),pp.maxwin,length(pp.chosen_g));
for pair=pp.chosen_p(pp.plotn_p)
    figure; subplot_dummy=0;
    overtitle=['Scatter Plot of ERPCOH for ',scl.p_label{pair}];
    for freq=pp.chosen_freq
        for win=1:pp.maxwin
            [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
            [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
            subplot_dummy=subplot_dummy+1;
            subplot(length(pp.chosen_freq),pp.maxwin,subplot_dummy);
            %
            for group=pp.chosen_g
            coh_scatterdata_x=squeeze(mean(mean(cohdata(t_start:t_end,freq,pp.cond_diff{1},pair,s_inds_g(:,group)),1),3));
            coh_scatterdata_y=squeeze(mean(mean(cohdata(t_start:t_end,freq,pp.cond_diff{2},pair,s_inds_g(:,group)),1),3));
            %test out a baseline subtraction
            %coh_scatterdata_x=squeeze(mean(mean(cohdata(t_start:t_end,freq,pp.cond_diff{1},pair,s_inds_g(:,scl.g_all)),1),3))-squeeze(mean(mean(cohdata(1:scl.t_zero,freq,pp.cond_diff{1},pair,s_inds_g(:,scl.g_all)),1),3));
            %coh_scatterdata_y=squeeze(mean(mean(cohdata(t_start:t_end,freq,pp.cond_diff{2},pair,s_inds_g(:,scl.g_all)),1),3))-squeeze(mean(mean(cohdata(1:scl.t_zero,freq,pp.cond_diff{1},pair,s_inds_g(:,scl.g_all)),1),3));
            %
            %note the subject indices of high coh values
            %s_high_coh(pair,freq,win,1:length(find(coh_scatterdata_x>.9 & coh_scatterdata_y>.9)))=find(coh_scatterdata_x>.9 & coh_scatterdata_y>.9);
            %
            p_tfwin(:,pair,freq,win,group)=polyfit(coh_scatterdata_x,coh_scatterdata_y,1);
            %
            scatter(coh_scatterdata_x,coh_scatterdata_y,scl.g_color{group}); hold on;
            plot(linspace(scatter_limits(1),scatter_limits(2),100),linspace(scatter_limits(1),scatter_limits(2),100)*p_tfwin(1,pair,freq,win,group)+p_tfwin(2,pair,freq,win,group),scl.g_color{group});
            end
            plot(linspace(scatter_limits(1),scatter_limits(2),100),linspace(scatter_limits(1),scatter_limits(2),100),'k--'); hold off;
            %title(sprintf('%d - %d, %d Hz, %s',pp.t_start_ms(win),pp.t_end_ms(win),round(scl.freqs(freq)),scl.p_label{pair}));
            %xlabel([scl.cond_label{pp.cond_diff{1}}])
            %ylabel([scl.cond_label{pp.cond_diff{2}}])
            axis([0 1 0 1])
        end
    end
    adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
    tightfig; dragzoom;
end
%unique(s_high_coh)
clear_plotassistvars

%% HISTOGRAM ERPCOH
% PAIRS

%seed_pairs=[13 14 16:23 25:30;31:34,40:41,44:45,48:49,55:60;61:76];
seed_pairs=[1:30;31:60;61:90];
seed_label={'FZ','CZ','PZ'};

for seed=1:3
figure; subplot_dummy=0;
for freq=pp.chosen_freq
    for win=1:pp.maxwin
        subplot_dummy=subplot_dummy+1;
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        %
        %coh_histdata=squeeze(mean(mean(cohdata(t_start:t_end,freq,pp.cond_diff{1},pair,s_inds_g(:,scl.g_all)),1),3) - ...
        %    mean(mean(cohdata(t_start:t_end,freq,pp.cond_diff{2},pair,s_inds_g(:,scl.g_all)),1),3));
        coh_histdata=squeeze(mean(mean(mean(cohdata(t_start:t_end,freq,pp.cond_diff{1},seed_pairs(seed,:),s_inds_g(:,scl.g_all)),1),3),4) - ...
            mean(mean(mean(cohdata(t_start:t_end,freq,pp.cond_diff{2},seed_pairs(seed,:),s_inds_g(:,scl.g_all)),1),3),4));
        %coh_histdata=squeeze(...
        %    (mean(mean(cohdata(t_start:t_end,freq,pp.cond_diff{1},pair,s_inds_g(:,scl.g_all)),1),3) - mean(mean(cohdata(1:scl.t_zero,freq,pp.cond_diff{1},pair,s_inds_g(:,scl.g_all)),1),3)) - ...
        %    (mean(mean(cohdata(t_start:t_end,freq,pp.cond_diff{2},pair,s_inds_g(:,scl.g_all)),1),3) - mean(mean(cohdata(1:scl.t_zero,freq,pp.cond_diff{2},pair,s_inds_g(:,scl.g_all)),1),3)));            
        coh_histdata_se=std(coh_histdata);
        subplot(length(pp.chosen_freq),pp.maxwin,subplot_dummy);
        hist(coh_histdata,pp.hist_nbin); axis(pp.hist_ax); hold on;
        h=findobj(gca,'Type','patch'); set(h,'FaceColor',[.7 .7 .7],'EdgeColor','w');
        plot(ones(1,100)*mean(coh_histdata),linspace(pp.hist_ax(3),pp.hist_ax(4),100),'k--'); hold on;
        plot(linspace(mean(coh_histdata)-std(coh_histdata),mean(coh_histdata)+std(coh_histdata),100), ...
            ones(1,100)*pp.hist_ax(4)/20,'k','LineWidth',2); hold off;
        %[f,x]=ksdensity(coh_histdata); plot(x,f);
        [~,p,~,stat]=ttest(zeros(sum(s_inds_g(:,scl.g_all)),1),coh_histdata);
        z=mean(coh_histdata)/coh_histdata_se;
        text(pp.p_loc(1),pp.p_loc(2),sprintf('z=%3.1f, t=%3.1f, p=%3.3f',z,stat.tstat,p));
        title(sprintf('%d - %d, %d Hz, %s',pp.t_start_ms(win),pp.t_end_ms(win),round(scl.freqs(freq)),seed_label{seed}));
        xlabel([[scl.cond_label{pp.cond_diff{1}}],' - ',[scl.cond_label{pp.cond_diff{2}}]])
    end
end
tightfig; dragzoom;
end
clear_plotassistvars

%% look at scl.freqs for all conditions averaged across seeds, averaged across all S's

%seed_pairs=[1:30;31:60;61:90];
%seed_pairs=[24;37;76];
%seed_pairs=[9 10 13 14 16:30;1:19;1:19];
%seed_pairs=[13 14 16:23 25:30;31:34,40:41,44:45,48:49,55:60;61:76];
%seed_label={'FZ','CZ','PZ'};
seed_pairs=[1:30;31:60;61:90;91:120];
seed_label={'F4','F3','P3','P4'};

for group=pp.chosen_g(pp.plotn_g)
for seed=1:4
    figure;
    for cond=pp.plotn_cond
        sp(cond)=subplot(sp_d(1),sp_d(2),cond);
        for freq=pp.chosen_freq
            if cond==imp.maxconds+1
                coh_plot_data = mean(mean(mean(cohdata(:,freq,pp.cond_diff{1},seed_pairs(seed,:),s_inds_g(:,group)),3),4),5) - ...
        mean(mean(mean(cohdata(:,freq,pp.cond_diff{2},seed_pairs(seed,:),s_inds_g(:,group)),3),4),5);
            else
                coh_plot_data=mean(mean(cohdata(:,freq,cond,seed_pairs(seed,:),s_inds_g(:,group)),4),5);
            end
            plot(coh_plot_data,'Color',scl.f_color(freq,:)); hold on; %,'Color',scl.s_color(pp.chosen_s,:)); hold on
        end
        if cond==imp.maxconds+1
            grid on; axis([scl.t_start scl.t_end pp.chosen_coh_limit_diff(1) pp.chosen_coh_limit_diff(2)]);
        else
            grid on; axis([scl.t_start scl.t_end pp.chosen_coh_limit(1) pp.chosen_coh_limit(2)]);
        end
        plot(ones(10,1)*scl.t_zero,linspace(0,1,10),'k--'); hold off;
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
        title([scl.g_label{group},'/',seed_label{seed},'/',scl.cond_label{cond}])
        clickableLegend(scl.f_nlabel(pp.chosen_freq))
    end
%linkaxes(sp(1:end-1))
tightfig;
end
end
%dragzoom; distFig('s','ext');
clear_plotassistvars

%% image coherence in time-freq at a chosen pair

sp_rowlabel={''};
sp_columnlabel=scl.cond_label;
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';
subplot_dims=sp_d;
%pp.chosen_p=1:14;

prestim_interestpairs=[92 96 98 100 104]; %102 108
%prestim_interestpairs=[112 119]; %110 115 117

pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),imp.maxpairs,length(pp.plotn_cond),2);
for group=pp.chosen_g(pp.plotn_g)
for pair=[7 8 36 37] %7 8 12 36 37 41] %prestim_interestpairs %pp.chosen_p(pp.plotn_p)
    pp.figdum=pp.figdum+1;
    figure(pp.figdum); subplot_dummy=0; 
    overtitle{pp.figdum}=sprintf('%s / %s',scl.p_label{pair},scl.g_label{group});
    %overtitle=sprintf('Phase Coherence between %s for %s',scl.p_label{pair},scl.g_label{group});
    for cond=pp.plotn_cond
        subplot_dummy=subplot_dummy+1;
        subplot(sp_d(1),sp_d(2),subplot_dummy)
        if cond==imp.maxconds+1
            %coh_tf_data=mean(atanh(cohdata(:,:,pp.cond_diff{1},pair,s_inds_g(:,group))),5) - ...
            %    mean(atanh(cohdata(:,:,pp.cond_diff{2},pair,s_inds_g(:,group))),5);
            coh_tf_data=mean(cohdata(:,:,pp.cond_diff{1},pair,s_inds_g(:,group)),5) - ...
                mean(cohdata(:,:,pp.cond_diff{2},pair,s_inds_g(:,group)),5);
        else
            %coh_tf_data=mean(atanh(cohdata(:,:,cond,pair,s_inds_g(:,group))),5);
            coh_tf_data=mean(cohdata(:,:,cond,pair,s_inds_g(:,group)),5);
        end
        contourf(fliplr(coh_tf_data)',pp.n_contour);
        shading flat
        %imagesc(mean(cohstats(:,:,cond,pair,s_inds_g(:,group)),5)');
        %axis([scl.t_start scl.t_end 3 imp.maxfreqs-2]);
        axis([scl.t_start scl.t_end 1 20]);
        v(group, pair,subplot_dummy,:) = caxis;
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
        set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
        %title([scl.p_label{pair},' / ',scl.cond_label{cond},' / ',scl.g_label{group}]); grid on;
        grid on
        hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
    end
    adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
    %adorn_plots({'conds','Loss-Gain'},{'Gain','Loss'},'Time (ms)','Frequency (Hz)',overtitle);
end
end
v(v==0)=NaN;
c(1)=nanmin(nanmin(nanmin(v(:,:,1:end-1,1)))); c(2)=nanmax(nanmax(nanmax(v(:,:,1:end-1,2))));
c_diff(1)=nanmin(nanmin(nanmin(v(:,:,end,1)))); c_diff(2)=nanmax(nanmax(nanmax(v(:,:,end,2))));
cmap=makecmap(c);
cmap_diff=makecmap(c_diff);
for fig=pp.figdum_init+1:pp.figdum
    figure(fig)
    for splot=1:subplot_dummy
        if splot==subplot_dummy
            subplot(sp_d(1),sp_d(2),splot); caxis([c_diff(1) c_diff(2)]);
            colormap(cmap_diff);
        else
            subplot(sp_d(1),sp_d(2),splot); caxis([c(1) c(2)]);
            colormap(cmap);
        end
        colorbar;
        freezeColors;
        cbfreeze;
    end
    plottitle(overtitle{fig});
    tightfig;
    set_print_size(20,8);
end
clear_plotassistvars

%% SCATTER ERPCOH WITH AGE

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='Age (yrs)';
y_plotlabel='ERPCOH';

itc_axes=[min(s_demogs.age_eeg) max(s_demogs.age_eeg) 0 0.8];
itc_axes_diff=[min(s_demogs.age_eeg) max(s_demogs.age_eeg) -0.5 0.5];

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_start_hz),pp.maxwin);
for pair=[1:12:120] %pp.chosen_p(pp.plotn_p)
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    figure; subplot_dummy=0;
    overtitle=sprintf('ERPCOH over age at %s in %1.1f - %1.1f Hz', ...
        scl.p_label{pair},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
    for cond=pp.plotn_cond
    for win=1:pp.maxwin
    [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
    %    
    for group=pp.chosen_g(pp.plotn_g)
    erpcoh_scatterdata_x=s_demogs.age_eeg(s_inds_g(:,group));
    if cond==imp.maxconds+1
        erpcoh_scatterdata_y=squeeze(mean(mean(mean(cohdata(...
            t_start:t_end,f_end:f_start,pp.cond_diff{1},pair,s_inds_g(:,group)),1),2),3)) -...
            squeeze(mean(mean(mean(cohdata(...
            t_start:t_end,f_end:f_start,pp.cond_diff{2},pair,s_inds_g(:,group)),1),2),3));
        axis(itc_axes_diff);
    else
        erpcoh_scatterdata_y=squeeze(mean(mean(mean(cohdata(...
            t_start:t_end,f_end:f_start,cond,pair,s_inds_g(:,group)),1),2),3));
        axis(itc_axes);
    end

    %p_tfwin(:,chan,freq_range,win)=robustfit(itc_scatterdata_x,itc_scatterdata_y);
    %
    scatter_h(group)=scatter(erpcoh_scatterdata_x,erpcoh_scatterdata_y, scl.g_color{group}); hold on;
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

%% ERPCOH - create condition bar plots with error bars

sp_rowlabel=make_freqlabels(pp.f_start_hz(pp.plotn_f),pp.f_end_hz(pp.plotn_f));
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='Group Condition-Means';
y_plotlabel='Phase Coherence';

width=[];
bw_xlabel=[];
bw_ylabel=[];
bw_colormap=bone; %pmkmp(length(barvalues),'cubicl');
gridstatus='y';
error_sides=1;
legend_type='plot';
bar_glabel={scl.g_label{pp.chosen_g(pp.plotn_g)}};
%bar_condlabel={'Go','NoGo'};
bar_condlabel={scl.cond_label{pp.plotn_cond(1:end-1)}};

%figures are chans, columns are time windows, rows are frequency bands
for pair=pp.chosen_p(pp.plotn_p)
figure;
subplot_dummy=0;
overtitle=scl.p_label{pair};
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    for win=1:length(pp.t_start_ms)
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        subplot_dummy=subplot_dummy+1;
        subplot(length(pp.plotn_f),pp.maxwin,subplot_dummy); %length(pp.f_start_hz),pp.maxwin,subplot_dummy);
        %
        gdum=0;
        for group=pp.chosen_g(pp.plotn_g)
        gdum=gdum+1;
        ranova_data{gdum}=squeeze(mean(mean(cohdata(t_start:t_end,f_end:f_start,pp.plotn_cond(1:end-1),pair,s_inds_g(:,group)),1),2))';
        bar_data(group,:)=squeeze(mean(mean(mean(cohdata(t_start:t_end,f_end:f_start,pp.plotn_cond(1:end-1),pair,s_inds_g(:,group)),1),2),5));
        %bar_data(group,:)=squeeze(mean(mean(mean(cohdata(t_start:t_end,f_end:f_start,:,pair,s_inds_g(:,group)),1),2),5) -...
        %    mean(mean(mean(cohdata(1:scl.t_zero,f_end:f_start,:,pair,s_inds_g(:,group)),1),2),5));
        bar_data_se(group,:)=std(squeeze(mean(mean(cohdata(t_start:t_end,f_end:f_start,pp.plotn_cond(1:end-1),pair,s_inds_g(:,group)),1),2)),0,2)/sqrt(sum(s_inds_g(:,group)));
        %bar_data_se(group,:)=std(squeeze(mean(mean(cohdata(t_start:t_end,f_end:f_start,:,pair,s_inds_g(:,group)),1),2) -...
        %    mean(mean(cohdata(1:scl.t_zero,f_end:f_start,:,pair,s_inds_g(:,group)),1),2)),0,2)/sqrt(sum(s_inds_g(:,group)));
        %
        end
        %bar_title=sprintf('%s, %1.1f - %1.1f, %d - %d ms',scl.p_label{pair},...
        %    pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),pp.t_start_ms(win),pp.t_end_ms(win));
        bar_title=[];
        bar_h=barweb(bar_data(pp.chosen_g,:),bar_data_se(pp.chosen_g,:),width,bar_glabel,...
            bar_title,bw_xlabel,bw_ylabel,bw_colormap,gridstatus,bar_condlabel,error_sides,legend_type);
        axis 'auto y'
        axis([0 4 0.2 .6])
        [p,ranova_table]=anova_rm(ranova_data,'off');
        text(2.5,0.35,sprintf('main p=%1.3f',p(1)))
        text(2.5,0.28,sprintf('group p=%1.3f',p(2)))
        %text(2.5,0.21,sprintf('int p=%1.3f',p(3)))
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,[length(pp.plotn_f),pp.maxwin]);
end
%linkaxes
%tightfig;
clear_plotassistvars

%% image the coherence statistic in time-freq at a chosen pair

%
pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),imp.maxpairs,length(pp.plotn_cond),2);
for group=pp.chosen_g(pp.plotn_g)
for pair=pp.chosen_p(pp.plotn_p)
    pp.figdum=pp.figdum+1;
    figure(pp.figdum); subplot_dummy=0; 
    for cond=pp.plotn_cond
        subplot_dummy=subplot_dummy+1;
        subplot(sp_d(1),sp_d(2),subplot_dummy)
        if cond==imp.maxconds+1
            coh_tf_data=mean(cohstats(:,:,pp.cond_diff{1},pair,s_inds_g(:,group)),5) - ...
                mean(cohstats(:,:,pp.cond_diff{2},pair,s_inds_g(:,group)),5);
        else
            coh_tf_data=mean(cohstats(:,:,cond,pair,s_inds_g(:,group)),5);
        end
        contourf(fliplr(coh_tf_data)');
        colormap(pmkmp(256,pp.pmkmp_scheme))
        shading flat
        %imagesc(mean(cohstats(:,:,cond,pair,s_inds_g(:,group)),5)');
        %axis([scl.t_start scl.t_end 3 imp.maxfreqs-2]);
        axis([scl.t_start scl.t_end 1 20]);
        v(group, pair,subplot_dummy,:) = caxis;
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); xlabel('Time (ms)');
        set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); ylabel('Frequency (Hz)');
        title([scl.p_label{pair},' / ',scl.cond_label{cond},'/',scl.g_label{group}]); grid on;
        hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
        colorbar;
    end
end
end
v(v==0)=NaN;
c(1)=nanmin(nanmin(nanmin(v(:,:,1:end-1,1)))); c(2)=nanmax(nanmax(nanmax(v(:,:,1:end-1,2))));
c_diff(1)=nanmin(nanmin(nanmin(v(:,:,end,1)))); c_diff(2)=nanmax(nanmax(nanmax(v(:,:,end,2))));
for fig=pp.figdum_init+1:pp.figdum
    figure(fig)
    for splot=pp.plotn_cond;
        if splot==imp.maxconds+1
            subplot(sp_d(1),sp_d(2),splot); caxis([c_diff(1) c_diff(2)]);
        else
            subplot(sp_d(1),sp_d(2),splot); caxis([c(1) c(2)]);
        end
    end
    tightfig;
end

%% image coherence in time-freq, averaged across seed pairs

%seed_pairs=[1:30;31:60;61:90];
seed_pairs=[1:30;31:60;61:90;91:120];
%seed_pairs=[1:30;31:60];
%seed_pairs=[13 14 16:23 25:30;31:34,40:41,44:45,48:49,55:60;61:76];
%seed_label={'C3','C4'};
seed_label={'F4','F3','P3','P4'};

%pp.cmap=colormap(hsv(64));
%pp.cmap=cbrewer('seq','BuGn',9);
%pp.cmap=pmkmp(256,pp.pmkmp_scheme);

pp.figdum=pp.figdum_init;

v=zeros(length(pp.chosen_g),length(seed_label),length(pp.plotn_cond),2);
for group=pp.chosen_g(pp.plotn_g)
for seed=1:length(seed_label)
    pp.figdum=pp.figdum+1;
    figure(pp.figdum); subplot_dummy=0; 
    for cond=pp.plotn_cond
        subplot_dummy=subplot_dummy+1;
        subplot(sp_d(1),sp_d(2),subplot_dummy)
        %
        if cond==imp.maxconds+1
            coh_seed_avg=mean(mean(mean(cohdata(:,:,pp.cond_diff{1},seed_pairs(seed,:),s_inds_g(:,group)),3),4),5)-...
                mean(mean(mean(cohdata(:,:,pp.cond_diff{2},seed_pairs(seed,:),s_inds_g(:,group)),3),4),5);
            %coh_seed_avg=(mean(mean(mean(cohdata(:,:,pp.cond_diff{1},seed_pairs(seed,:),s_inds_g(:,group)),3),4),5) - ...
            %    repmat(squeeze(mean(mean(mean(mean(cohdata(1:scl.t_zero,:,pp.cond_diff{1},seed_pairs(seed,:),s_inds_g(:,group)),1),3),4),5)),imp.maxtimepts,1)) -...
            %    (mean(mean(mean(cohdata(:,:,pp.cond_diff{2},seed_pairs(seed,:),s_inds_g(:,group)),3),4),5) - ...
            %    repmat(squeeze(mean(mean(mean(mean(cohdata(1:scl.t_zero,:,pp.cond_diff{2},seed_pairs(seed,:),s_inds_g(:,group)),1),3),4),5)),imp.maxtimepts,1));
        else
            coh_seed_avg=mean(mean(cohdata(:,:,cond,seed_pairs(seed,:),s_inds_g(:,group)),4),5);
            %coh_seed_avg=mean(mean(cohdata(:,:,cond,seed_pairs(seed,:),s_inds_g(:,group)),4),5) - ...
                repmat(squeeze(mean(mean(mean(cohdata(1:scl.t_zero,:,cond,seed_pairs(seed,:),s_inds_g(:,group)),1),4),5)),imp.maxtimepts,1);
        end
        %
        contourf(fliplr(coh_seed_avg)',pp.n_contour);
        shading flat
        %
        %imagesc(coh_seed_avg');
        %
        colormap(pp.cmap)
        %
        axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
        v(group,seed,cond,:) = caxis;
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); xlabel('Time (ms)');
        set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); ylabel('Frequency (Hz)');
        title([seed_label{seed},' / ',scl.cond_label{cond},' / ',scl.g_label{group}])
        hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
    end
end
end
v(v==0)=NaN;
c(1)=nanmin(nanmin(nanmin(v(:,:,1:end-1,1)))); c(2)=nanmax(nanmax(nanmax(v(:,:,end-1,2))));
diff_cmin=nanmin(nanmin(v(:,:,end,1))); diff_cmax=nanmax(nanmax(v(:,:,end,2)));
for fig=pp.figdum_init+1:pp.figdum
    figure(fig)
    for splot=pp.plotn_cond
        if splot==imp.maxconds+1
            subplot(sp_d(1),sp_d(2),splot); caxis([diff_cmin diff_cmax]);
        else
            subplot(sp_d(1),sp_d(2),splot); caxis([c(1) c(2)]);
        end
        colorbar;
    end
    tightfig;
end

%% plot coherence pairs as lines on a topoplot
data=zeros(imp.maxchans,1);
h=figure;
topoplot(data,chan_locs,'style','blank');
hold on
for pair=1:imp.maxpairs
    line([chan_locs(param_struct.coherence_pairs(pair,1)).topo_x chan_locs(param_struct.coherence_pairs(pair,2)).topo_x],...
        [chan_locs(param_struct.coherence_pairs(pair,1)).topo_y chan_locs(param_struct.coherence_pairs(pair,2)).topo_y],...
        [2.1-randn*.05 2.1+randn*.05],'Color',scl.p_color(pair,:)); hold off; %disp(num2str(pair))
end

%% coherence in frequency band as topoplot with colored lines indicating strength, multiple time windows
sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='Time Windows';
y_plotlabel='Conditions';

%define re-scaling constants
linescale=[1,256];

%choose pair sub-set
cbar_ticks=5;

%line_limit=[-1 1];
%line_limit_diff=[-1 1];
%
line_limit=[0 .12];
line_limit_diff=[-.03 .06];

cmap=makecmap(line_limit);
cmap_diff=makecmap(line_limit_diff);

alpha=.000001;
plot_sig=true;

coh_linescale_mat=zeros(length(pp.plotn_cond),length(pp.chosen_g),length(pp.t_start_ms),length(pp.f_start_hz),imp.maxpairs);
dummydata=ones(length(chan_locs),1)*0.3;
pp.figdum=pp.figdum_init;
for group=pp.chosen_g(pp.plotn_g)
for freq_range=6 %pp.plotn_f
    %convert scl.freqs to points
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    pp.figdum=pp.figdum+1;
    figure(pp.figdum); subplot_dummy=0;
    overtitle{pp.figdum}=sprintf('Coherence Pairs, %s / %1.1f - %1.1f Hz',scl.g_label{group},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        subplot_dummy=subplot_dummy+1;
        subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy)
        if win~=pp.maxwin
        topoplot(dummydata,chan_locs,'style','blank','maplimits',[0 1],'electrodes','off'); hold on;
        for pair=1:imp.maxpairs
            %determine significance before further plotting
            if plot_sig
            gdum=0;
            for statsgroup=pp.chosen_g(pp.plotn_g)
                gdum=gdum+1;
                ranova_data{gdum}=squeeze(mean(mean(cohdata(t_start:t_end,f_end:f_start,:,pair,s_inds_g(:,statsgroup)),1),2))';
            end
            [p,~]=anova_rm(ranova_data,'off');
            %if ~any(p([1:2,4])<alpha)
            if p(1)>alpha
                continue
            end
            end
            %define [x1 x2], [y1 y2], and [z1 z2] of the arc
            x=[chan_locs(param_struct.coherence_pairs(pair,1)).topo_x chan_locs(param_struct.coherence_pairs(pair,2)).topo_x];
            y=[chan_locs(param_struct.coherence_pairs(pair,1)).topo_y chan_locs(param_struct.coherence_pairs(pair,2)).topo_y];
            %z=[2.1-randn*.05 2.1+randn*.05];
            %determine strength of coherence
            if cond==imp.maxconds+1
                paircoh_for_linecolor=mean(mean(mean(mean(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{1},pair,s_inds_g(:,group)),1),2),3),5) - ...
                    mean(mean(mean(mean(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{2},pair,s_inds_g(:,group)),1),2),3),5);
                coh_linescale_mat(cond,group,win,freq_range,pair)=paircoh_for_linecolor;
                paircoh_color=(paircoh_for_linecolor-line_limit_diff(1))/(line_limit_diff(2) - line_limit_diff(1))*(linescale(2)-1);
                paircoh_color = ceil(paircoh_color)+1;
                linecolor=cmap_diff(paircoh_color,:);
                linesize=abs(paircoh_for_linecolor)/range(line_limit_diff)*6;
            else
                %paircoh_for_linecolor=mean(mean(mean(cohdata(t_start:t_end,f_end:f_start,cond,pair,s_inds_g(:,group)),1),2),5);
                paircoh_for_linecolor=mean(mean(mean(cohdata(t_start:t_end,f_end:f_start,cond,pair,s_inds_g(:,group)),1),2),5) -...
                    mean(mean(mean(cohdata(1:scl.t_zero,f_end:f_start,cond,pair,s_inds_g(:,group)),1),2),5); %subtract baseline
                %paircoh_for_linecolor=mean(mean(mean(cohdata(t_start:t_end,f_end:f_start,cond,pair,s_inds_g(:,group)),1),2),5) -...
                %    mean(mean(mean(cohdata(scl.t_zero+1:imp.maxtimepts,f_end:f_start,cond,pair,s_inds_g(:,group)),1),2),5);
                %    %subtract post-stim coherence??
                coh_linescale_mat(cond,group,win,freq_range,pair)=paircoh_for_linecolor;
                paircoh_color=(paircoh_for_linecolor-line_limit(1))/(line_limit(2) - line_limit(1))*(linescale(2)-1);
                paircoh_color = ceil(paircoh_color)+1;
                linecolor=cmap(paircoh_color,:);
                linesize=abs(paircoh_for_linecolor)/range(line_limit)*6;
            end
            %store color for scaling
            %scale post hoc
            %paircoh_color = ceil ( ( paircoh_color - linescale(1) + 1 ) * linescale(2) / ( linescale(2) - linescale(1) + 1 ) );
            
            %define its color and width based on its strength
            %linecolor=pp.cmap_line(paircoh_color,:);
            %linesize=(paircoh_color/256)*5;
            %direction=mod(pair,2);
            direction=1;
            %plot the arc
            Draw_Arc_Clockwise([x(1),y(1)], [x(2),y(2)], linecolor, linesize, direction);
            hold on;
        end
        else
            set(gca,'Visible','Off');
        end
        hold off; %title(sprintf('%s / %s, %d - %d ms, %1.1f - %1.1f Hz',scl.g_label{group},scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),pp.f_start_hz(freq_range),pp.f_end_hz(freq_range)))
    end
end
end
end
for fig=pp.figdum_init+1:pp.figdum
figure(fig);
for cond=1:length(pp.plotn_cond)
    subplot(length(pp.plotn_cond),length(pp.t_start_ms),cond*length(pp.t_start_ms))
    cb_pos = get(gca,'Position') + [0.08 0 0 0];
    cb_pos(3) = 0.05; %set width
    if cond==imp.maxconds+1
        colormap(cmap_diff);
        h=colorscale([1 256], line_limit_diff, range(line_limit_diff)/5, 'vert','Position',cb_pos);
        %cb_ticks=cellstr(num2str(linspace(line_limit_diff(1),line_limit_diff(2),cbar_ticks+1)','%1.2f'))';
        %set(h,'YTickLabel',{cb_ticks{2:end}});
        freezeColors;
    else
        colormap(cmap);
        cb_lim = line_limit;
        h=colorscale([1 256], line_limit, range(line_limit)/5, 'vert','Position',cb_pos);
        %cb_ticks=cellstr(num2str(linspace(line_limit(1),line_limit(2),cbar_ticks+1)','%1.2f'))';
        %set(h,'YTickLabel',{cb_ticks{2:end}});
        freezeColors;
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle{fig},[length(pp.plotn_cond),length(pp.t_start_ms)]);
tightfig;
end
coh_linescale_mat(coh_linescale_mat==0)=NaN;
c=[nanmin(nanmin(nanmin(nanmin(nanmin(coh_linescale_mat(1:end-1,:,:,:,:)))))) nanmax(nanmax(nanmax(nanmax(nanmax(coh_linescale_mat(1:end-1,:,:,:,:))))))];
c_diff=[nanmin(nanmin(nanmin(nanmin(nanmin(coh_linescale_mat(end,:,:,:,:)))))) nanmax(nanmax(nanmax(nanmax(nanmax(coh_linescale_mat(end,:,:,:,:))))))];
clear_plotassistvars

%% image coherence as a topoplot based on seeds

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

sp_rowlabel=scl.cond_label;
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel='Time Windows';
y_plotlabel='Conditions';

%choose pair sub-set
%seed_pairs=[7 13 14 16:23 25:30;31:34,40:41,44:45,48:49,55:60;61:76];
seed_pairs=[1:30;31:60;61:90;91:120];
%seed_pairs=[1:60;61:120];
%seed_label={'FPZ','OZ'};
seed_label={'F4','F3','P3','P4'};
seed_inds=[8,9,23,24]; %38, 58]; %8,9,23,24]; %17 18
maxseeds=length(seed_inds);

coh_topo_scale=[-.01 .08];
coh_diff_limits=[-.03 .02];
cmap=makecmap(coh_topo_scale);
cmap_diff=makecmap(coh_diff_limits);

seed_chanlocs=load('/export/home/mike/matlab/origin/coords/31chans_ns.mat');
seed_chanlocs=getfield(seed_chanlocs,'chan_locs');

%coh_topo_scale=[-1 1];
%coh_diff_limits=[-1 1];

v=zeros(length(pp.plotn_cond),length(pp.plotn_f),length(pp.chosen_g),maxseeds,pp.maxwin,2);
%pp.figdum=pp.figdum_init;
pp.figdum_init=pp.figdum;

for freq_range=1 %pp.plotn_f

%convert scl.freqs to points
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));

%figure('units','normalized','position',[.1 .1 .9 .4]);
for group=pp.chosen_g
for seed=[1 2 3 4]
pp.figdum=pp.figdum+1;
subplot_dummy=0; fig(pp.figdum);
overtitle{pp.figdum}=sprintf('Event-Related Coherence with %s, %s / %1.1f - %1.1f Hz',seed_label{seed},scl.g_label{group},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
for cond=pp.plotn_cond
    for win=1:pp.maxwin
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        subplot_dummy=subplot_dummy+1;
        subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy)
        % grab the relevant data
        if cond==imp.maxconds+1
            if sum(s_inds_g(:,group))<600
            coh_topo_data=meanx(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{1},seed_pairs(seed,:),s_inds_g(:,group)),4) -...
                meanx(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{2},seed_pairs(seed,:),s_inds_g(:,group)),4);
            else
            coh_topo_data=zeros(length(seed_pairs(seed,:)),1);
            sdum=0;
            for seedpair=seed_pairs(seed,:)
            sdum=sdum+1;
            coh_topo_data(sdum)=meanx(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{1},seedpair,s_inds_g(:,group)),4) -...
                meanx(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{2},seedpair,s_inds_g(:,group)),4);
            end
            end
            coh_topo_data=insertrows(coh_topo_data,0,seed_inds(seed)-1);
        else
            if sum(s_inds_g(:,group))<600
            coh_topo_data=meanx(cohdata(t_start:t_end,f_end:f_start,cond,seed_pairs(seed,:),s_inds_g(:,group)),4) - ...
                meanx(cohdata(1:scl.t_start,f_end:f_start,pp.cond_diff{1},seed_pairs(seed,:),s_inds_g(:,group)),4);
            else
            coh_topo_data=zeros(length(seed_pairs(seed,:)),1);
            sdum=0;
            for seedpair=seed_pairs(seed,:)
            sdum=sdum+1;
            coh_topo_data(sdum)=meanx(cohdata(t_start:t_end,f_end:f_start,cond,seedpair,s_inds_g(:,group)),4) - ...
                meanx(cohdata(1:scl.t_start,f_end:f_start,pp.cond_diff{1},seedpair,s_inds_g(:,group)),4);
            end
            end
            %coh_topo_data=zeros(30,1);
            coh_topo_data=insertrows(coh_topo_data,0,seed_inds(seed)-1);
        end
        v(cond,freq_range,group,seed,win,1)=min(coh_topo_data);
        v(cond,freq_range,group,seed,win,2)=max(coh_topo_data);
        if win~=pp.maxwin
        if cond==imp.maxconds+1
            topoplot(coh_topo_data,seed_chanlocs,'maplimits',[coh_diff_limits(1) coh_diff_limits(2)],'electrodes','off',...
                'colormap',cmap_diff,'style','fill','numcontour',7);
            freezeColors;
        else
            topoplot(coh_topo_data,seed_chanlocs,'maplimits',[coh_topo_scale(1) coh_topo_scale(2)],'electrodes','off',...
                'colormap',cmap,'style','fill','numcontour',7);
            freezeColors;
        end
        shading flat
        else
            set(gca,'Visible','Off');
        end
    end
    %color bar
    cb_pos = get(gca,'Position') + [0.08 0 0 0];
    cb_pos(3) = 0.05; %set width
    if cond==imp.maxconds+1
        cb_lim = coh_diff_limits;
    else
        cb_lim = coh_topo_scale;
    end
    colorscale([1 256], cb_lim, range(cb_lim)/5, 'vert', ...
        'Position',cb_pos);
    freezeColors;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel, ...
    overtitle{pp.figdum},[length(pp.plotn_cond),length(pp.t_start_ms)]);
fig(pp.figdum, 'units', 'centimeters', 'width', 17, 'height', length(pp.plotn_cond)*3, ...
    'font', 'Monospaced', 'fontsize', 11);
set_print_size(17,17/scl.phi);
end
end
end

c(1)=min(min(min(min(min(v(1:end-1,:,:,:,:,1)))))); c(2)=max(max(max(max(max(v(1:end-1,:,:,:,:,2))))));
c_diff(1)=min(min(min(min(v(end,:,:,:,:,1))))); c_diff(2)=max(max(max(max(v(end,:,:,:,:,2)))));
clear_plotassistvars

%% compare coherence and ITC

 %14 is FZ-CZ
chosen_chan=7;

for pair=pp.chosen_p(pp.plotn_p)
figure

for group=pp.chosen_g
chan_coh_measure=squeeze(mean(mean(cohdata(180:231,14:16,1,pair,s_inds_g(:,group)),1),2));
chan_coh_stat=squeeze(mean(mean(cohstats(180:231,14:16,1,pair,s_inds_g(:,group)),1),2));
trial_coh_measure=squeeze(mean(mean(itcdata(180:231,chosen_chan,14:16,1,s_inds_g(:,group)),1),3));

scatter3(trial_coh_measure,chan_coh_measure,chan_coh_stat,scl.g_color{group}); hold on;
end
grid on;
xlabel('ITC'); ylabel('Channel Cross-Coherence'); zlabel('Z Score');
title([scl.p_label{pair},' / ', scl.chan_label{chosen_chan}]);
end

%% scatter coherence in two TF-window/pair combos
% does pre-stimulus alpha coherence predict other measures in individuals?

xpair=[108 102];
xwin_t=[-200 0];
xwin_f=[9 11];

ypair=[8 8];
ywin_t=[200 400];
ywin_f=[4.6 6.4];

[~,xt_s]=min(abs(scl.t_ms-xwin_t(1))); [~,xt_e]=min(abs(scl.t_ms-xwin_t(2)));
[~,xf_s]=min(abs(scl.freqs-xwin_f(1))); [~,xf_e]=min(abs(scl.freqs-xwin_f(2)));

[~,yt_s]=min(abs(scl.t_ms-ywin_t(1))); [~,yt_e]=min(abs(scl.t_ms-ywin_t(2)));
[~,yf_s]=min(abs(scl.freqs-ywin_f(1))); [~,yf_e]=min(abs(scl.freqs-ywin_f(2)));

figure; subplot_dummy=0;
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(sp_d(1),sp_d(2),subplot_dummy)
    gdum=0;
    
for group=pp.chosen_g
    
    gdum=gdum+1;
    if cond==imp.maxconds+1
    xdata=meanx(cohdata(xt_s:xt_e,xf_e:xf_s,pp.cond_diff{1},xpair(gdum),s_inds_g(:,group)),5) -...
        meanx(cohdata(xt_s:xt_e,xf_e:xf_s,pp.cond_diff{2},xpair(gdum),s_inds_g(:,group)),5);
    ydata=meanx(cohdata(yt_s:yt_e,yf_e:yf_s,pp.cond_diff{1},ypair(gdum),s_inds_g(:,group)),5) - ...
        meanx(cohdata(yt_s:yt_e,yf_e:yf_s,pp.cond_diff{2},ypair(gdum),s_inds_g(:,group)),5);
    else
    xdata=meanx(cohdata(xt_s:xt_e,xf_e:xf_s,cond,xpair(gdum),s_inds_g(:,group)),5);
    ydata=meanx(cohdata(yt_s:yt_e,yf_e:yf_s,cond,ypair(gdum),s_inds_g(:,group)),5);
    end
    
    scatter(xdata,ydata,scl.g_color{group}); hold on;
end
hold off;
if cond==imp.maxconds+1
axis([-.5 .5 -.5 .5]); grid on;
else
axis([0 1 0 1]); grid on;
end
end

%% plot per-subject lines of the shape of results for a certain measure across conditions (ERPCOH)

%chosen_pair=[5 7 8 12 36 37 41]
chosen_pair=8;
win_t=[200 400];
%win_t=[200 300];
win_f=[4 5.3]; %[4 5.3];

colorkey={'r','g','b','m','k','c'};
sitekey={'UConn','Indiana','Iowa','SUNY','WashU','UCSD'};
counter=zeros(1,6);

[~,t_s]=min(abs(scl.t_ms-win_t(1))); [~,t_e]=min(abs(scl.t_ms-win_t(2)));
[~,f_s]=min(abs(scl.freqs-win_f(1))); [~,f_e]=min(abs(scl.freqs-win_f(2)));

figure; gdum=0;
overtitle=sprintf('ERPCOH at %s in %1.1f - %1.1f Hz from %d - %d ms',...
    scl.p_label{chosen_pair},win_f(1),win_f(2),win_t(1),win_t(2));
for group=pp.chosen_g(pp.plotn_g)
    conddiff_dir=zeros(imp.maxconds,1);
    gdum=gdum+1;
    subplot(1,2,gdum)
    linedata=meanx(cohdata(t_s:t_e,f_e:f_s,:,chosen_pair,s_inds_g(:,group)),[3 5]);
    %linedata=meanx(cohdata(t_s:t_e,f_e:f_s,:,chosen_pair,s_inds_g(:,group)),[3 5]) - ...
    %    meanx(cohdata(1:scl.t_start,f_e:f_s,:,chosen_pair,s_inds_g(:,group)),[3 5]);
    ranova_data{gdum}=linedata';
    for s=1:size(linedata,2)
    %site=str2double(mat_list{s}(65));
    %counter(site)=counter(site)+1;
    if linedata(pp.cond_diff{1},s) > linedata(pp.cond_diff{2},s)
        mark='b-o'; conddiff_dir(1)=conddiff_dir(1)+1;
    else
        mark='m--*'; conddiff_dir(2)=conddiff_dir(2)+1;
    end
    %subplot(1,6,site);
    %title(sitekey{site});
    plot(linedata(:,s),mark); hold on; axis([0 imp.maxconds+1 -.2 1]);
    end
    conddiff_text=sprintf(['%s > %s : %d/%d (%2.0f%%), M=%0.2f',char(177),'%0.2f'],...
        scl.cond_label{pp.cond_diff{1}},scl.cond_label{pp.cond_diff{2}},...
        conddiff_dir(1),sum(s_inds_g(:,group)),...
        conddiff_dir(1)/sum(s_inds_g(:,group))*100,...
        abs(diff(mean(linedata(:,linedata(pp.cond_diff{1},:)>...
        linedata(pp.cond_diff{2},:)),2))),...
        abs(diff(std(linedata(:,linedata(pp.cond_diff{1},:)>...
        linedata(pp.cond_diff{2},:)),0,2))));
    text(0.5,0.8,conddiff_text)
    set(gca,'XTick',[1:imp.maxconds],'XTickLabel',scl.cond_label);
    xlabel('Condition');
    ylabel('ERPCOH');
    title(scl.g_label{group});    
end
[p,ranova_table]=anova_rm(ranova_data,'off');
plottitle(overtitle);
clear_plotassistvars