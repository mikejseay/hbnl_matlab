%% sample histograms

figure;
for group=pp.chosen_g(pp.plotn_g)
    subplot(2,1,group-1)
    hist(s_demogs.age_eeg(s_inds_g(:,group)),18:3:45);
    axis([18 45 0 30]); grid on;
end
set_print_size(8,16);

%% check ERP peak data

figure; subplot_dummy=0;
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    sp(subplot_dummy)=subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy);
    for group=pp.chosen_g
        if cond==imp.maxconds+1
        scatter(peakmat(s_inds_g(:,group),pp.cond_diff{1},2) - peakmat(s_inds_g(:,group),pp.cond_diff{2},2),...
            peakmat(s_inds_g(:,group),pp.cond_diff{1},1) - peakmat(s_inds_g(:,group),pp.cond_diff{2},1),scl.g_color{group}); hold on;
        else
        scatter(peakmat(s_inds_g(:,group),cond,2),...
            peakmat(s_inds_g(:,group),cond,1),scl.g_color{group}); hold on;
        end
    end
    hold off; set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms);
    axis tight; title(sprintf('%s',scl.cond_label{cond}));
    xlabel('Peak Latency (ms)'); ylabel('Peak Amplitude (uV)');
end
%linkaxes(sp)
clear_plotassistvars

%% plot ERPS with groups superimposed

%warning off MATLAB:linkaxes:RequireDataAxes

pp.figdum=pp.figdum_init;

% [20 22 24 26 28 49 51 59] %
% cyan = G1, magenta = G2
lims = [-8 19];

h_line=zeros(length(pp.plotn_cond),length(pp.chosen_g));
for chan=[7 16 25] %pp.chosen_chan(pp.plotn_chan)
    pp.figdum=pp.figdum+1;
    figure(pp.figdum)
    subplot_dummy=0;
    for cond=pp.plotn_cond
        subplot_dummy=subplot_dummy+1;
        sp(cond)=subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy);
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
            h=shadedErrorBar(1:imp.erpmaxtimepts,erp_plot_data,erp_plot_data_std,scl.g_color{group}); hold on;
            h_line(cond,group)=h.mainLine;
        end
        axis([scl.t_start_ind_erp imp.erpmaxtimepts lims(1) lims(2)]);
        %axis tight
        vline(scl.t_zero_erp,'k--'); hold off;
        set(gca,'XTick',scl.t_xtick_erp,'XTickLabel',scl.t_xtick_ms)
        grid on
        title([scl.chan_label{chan},'/',scl.cond_label{cond}])
    end
tightfig;
%linkaxes(sp(1:end-1))
%set_print_size(17,17/scl.phi);
set(gcf,'position',[120 120 987 300]);
end
clear_plotassistvars

%% plot group ERPs from avg_h1 data structure

% pp.figdum=pp.figdum_init;

% [20 22 24 26 28 49 51 59] %
% cyan = G1, magenta = G2
lims = [-8 19];

h_line=zeros(length(pp.plotn_cond),length(pp.chosen_g));
for chan=[7 16 25] %pp.chosen_chan(pp.plotn_chan)
    pp.figdum=pp.figdum+1;
    figure(pp.figdum)
    subplot_dummy=0;
    for cond=pp_h1.plotn_cond
        subplot_dummy=subplot_dummy+1;
        sp(cond)=subplot(pp_h1.sp_d(1),pp_h1.sp_d(2),subplot_dummy);
        for group=pp.chosen_g(pp.plotn_g)
            if cond==imp.maxconds+1
                erp_plot_data = mean(mean(erpdata_h1(:,chan,pp.cond_diff{1},s_inds_g(:,group)),3),4) - ...
            mean(mean(erpdata_h1(:,chan,pp.cond_diff{2},s_inds_g(:,group)),3),4);
                erp_plot_data_std = std(mean(erpdata_h1(:,chan,pp.cond_diff{1},s_inds_g(:,group)),3) - ...
            mean(erpdata_h1(:,chan,pp.cond_diff{2},s_inds_g(:,group)),3),0,4)/sqrt(sum(s_inds_g(:,group)));
            else
                erp_plot_data=mean(erpdata_h1(:,chan,cond,s_inds_g(:,group)),4);
                erp_plot_data_std=std(erpdata_h1(:,chan,cond,s_inds_g(:,group)),0,4)/sqrt(sum(s_inds_g(:,group)));
            end
            %plot(ero_plot_data,scl.g_color{group}); hold on; %,'Color',scl.s_color(pp.chosen_s,:)); hold on
            h=shadedErrorBar(1:scl_h1.erpmaxtimepts,erp_plot_data,erp_plot_data_std,scl.g_color{group}); hold on;
            h_line(cond,group)=h.mainLine;
        end
        axis([scl_h1.t_start_ind_erp scl_h1.t_end_ind_erp lims(1) lims(2)]);
        %axis tight
        vline(scl_h1.t_zero_erp,'k--'); hold off;
        set(gca,'XTick',scl_h1.t_xtick_erp,'XTickLabel',scl.t_xtick_ms)
        grid on
        title([scl.chan_label{chan},'/',scl_h1.cond_label{cond}])
    end
tightfig;
%linkaxes(sp(1:end-1))
%set_print_size(17,17/scl.phi);
set(gcf,'position',[120 120 987 300]);
end
clear_plotassistvars


%% plot ERPs with conditions superimposed

cond_colors={'r','g','b','m','r--','g--','b--','m--','k--','c'};
for chan=[7 16 25] %pp.chosen_chan
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
        vline(scl.t_zero_erp,'k--'); hold off;
        set(gca,'XTick',scl.t_xtick_erp,'XTickLabel',scl.t_xtick_ms)
        title([scl.chan_label{chan},'/',scl.g_label{group}])
        clickableLegend(scl.cond_label(pp.plotn_cond))
    end
tightfig;
linkaxes(sp)
end
set_print_size(18,9/scl.phi)
clear_plotassistvars

%% plot ERPs with conditions superimposed SPECIAL ANNOTATED

t1 = -1500; %choice stim
t2 = -1000; %response
[~,t1_pt] = min(abs(t1 - scl.t_ms_erp));
[~,t2_pt] = min(abs(t2 -scl.t_ms_erp));

%cond_colors={'r','r--','g--','g','k'};
cond_colors={'r','g','k--'};
for chan=7 %pp.chosen_chan
    figure;
    hold on;
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
            plot(erp_plot_data,cond_colors{cond}); 
            % dashline
        end
        axis tight
        vline(scl.t_zero_erp,'k');
        vline(t1_pt,'k');
        vline(t2_pt,'k'); hold off;
        set(gca,'XTick',scl.t_xtick_erp(1:2:end),'XTickLabel',scl.t_xtick_ms(1:2:end))
        title([scl.chan_label{chan},'/',scl.g_label{group}])
        clickableLegend(scl.cond_label(pp.plotn_cond))
    end
    tightfig;
    linkaxes(sp)
    %set(gca,'Visible','Off');
end
set_print_size(18,9/scl.phi)
clear_plotassistvars

%% plot ERPs / topographies for individual subjects with mike cohen's ERPviewer

suspects = [5 6 20 23 29 33 39 56 60];

subject = 11;

for s=subject %1:60 %suspects
    %disp(s)
for cond=pp.plotn_cond
    if cond==imp.maxconds+1
        erpviewerx(scl.t_ms_erp(scl.t_start_ind_erp:imp.erpmaxtimepts), ...
            squeeze(erpdata(scl.t_start_ind_erp:imp.erpmaxtimepts,:,pp.cond_diff{1},s))' - ...
            squeeze(erpdata(scl.t_start_ind_erp:imp.erpmaxtimepts,:,pp.cond_diff{2},s))', ...
            chan_locs);
    else
        h(cond) = erpviewerx(scl.t_ms_erp(scl.t_start_ind_erp:imp.erpmaxtimepts), ...
            squeeze(erpdata(scl.t_start_ind_erp:imp.erpmaxtimepts,:,cond,s))', ...
            chan_locs, ...
            [scl.cond_label{cond},'-',num2str(n_trials_all(cond,s))]);
    end
end
%pause
%close all hidden
end

%%

for s=subject %1:60 %suspects
    %disp(s)
for cond=pp_h1.plotn_cond
    if cond==imp.maxconds+1
        erpviewerx(scl_h1.t_ms_erp(scl_h1.t_start_ind_erp:scl_h1.t_end_ind_erp), ...
            squeeze(erpdata_h1(scl_h1.t_start_ind_erp:scl_h1.t_end_ind_erp,:,pp_h1.cond_diff{1},s))' - ...
            squeeze(erpdata_h1(scl_h1.t_start_ind_erp:scl_h1.t_end_ind_erp,:,pp_h1.cond_diff{2},s))', ...
            chan_locs);
    else
        h(cond) = erpviewerx(scl_h1.t_ms_erp(scl_h1.t_start_ind_erp:scl_h1.t_end_ind_erp), ...
            squeeze(erpdata_h1(scl_h1.t_start_ind_erp:scl_h1.t_end_ind_erp,:,cond,s))', ...
            chan_locs, ...
            [scl_h1.cond_label{cond},'-']);
    end
end
%pause
%close all hidden
end

%% do your own butterfly plot / global field power plot


lims = [-15 30];

%s_inds=26;
%s_inds=s_inds_g(1,:);
s_inds=find(s_inds_g(:,1))';

figure;
for s=s_inds(1:end)
clear yn
clf
subplot_dummy=0;
for cond=pp.plotn_cond
    %butterfly
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.plotn_cond),2,subplot_dummy);
    if cond==imp.maxconds+1
        erp_plot_data = meanx(erpdata(:,:,pp.cond_diff{1},s),[1 2]) - ...
            meanx(erpdata(:,:,pp.cond_diff{2},s),[1 2]);
    else
        erp_plot_data = meanx(erpdata(:,:,cond,s),[1 2]);
    end
    plot(erp_plot_data); %times are row, channels columns (check)
    axis([scl.t_start_erp scl.t_end_erp lims(1) lims(2)]);
    vline(scl.t_zero_erp,'k--'); hold off;
    set(gca,'XTick',scl.t_xtick_erp,'XTickLabel',scl.t_xtick_ms)

    %topographical variance
    subplot_dummy=subplot_dummy+1;
    sp(subplot_dummy)=subplot(length(pp.plotn_cond),2,subplot_dummy);
    if cond==imp.maxconds+1
        gfp_plot_data = var( meanx(erpdata(:,:,pp.cond_diff{1},s),[1 2]) - ...
            meanx(erpdata(:,:,pp.cond_diff{2},s),[1 2]) , 0, 2 );
        %gfp_plot_data = var( meanx(erpdata(:,:,pp.cond_diff{1},:),[1 2]),0,2) - ...
        %    var(meanx(erpdata(:,:,pp.cond_diff{2},:),[1 2]),0,2);
    else
        gfp_plot_data = var( meanx(erpdata(:,:,cond,s),[1 2]), 0, 2 );
    end
    plot(gfp_plot_data);
    axis([scl.t_start_erp scl.t_end_erp 0 100]);
    vline(scl.t_zero_erp,'k--'); hold off;
    set(gca,'XTick',scl.t_xtick_erp,'XTickLabel',scl.t_xtick_ms)
end
%linkaxes(sp);
plottitle(num2str(s));

yn=input(sprintf('S%d:',s),'s');
if strcmpi(yn,'n')
    rejmat(s)=true;
else
    rejmat(s)=false;
end
    
end

%% all encompassing single-subject data check figure
% butterfly plot for conditions and difference with major electrodes bolded
% GFP for conditions
% topographies at key time-point (~250 ms) for conditions and differences
% ERSP / ITC at key channel (Fz)
% ERCOH for key relation (intra-frontal)

s_inds=find(s_inds_g(:,3))';
%rejmat=false(size(s_inds));

%key time window
ck_win=[230 330];
[~,t_start_ck]=min(abs(scl.t_ms_erp-ck_win(1)));
[~,t_end_ck]=min(abs(scl.t_ms_erp-ck_win(2)));

%key channel
ck_chan=7;

%key regional hypothesis
ck_region=1;
plot_hypinds=find(opt.pair_inds==ck_region);

%baseline
[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

cmap=makecmap([-1 1],0);

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel={'Butterfly','GFP','Topo in ck_win','ERSP','ITPC','IRPS'};
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),6];

figure;
for s=s_inds
clear yn
clf
subplot_dummy=0;
for cond=pp.plotn_cond
    if cond==imp.maxconds+1
        erp_plot_data = meanx(erpdata(:,:,pp.cond_diff{1},s),[1 2]) - ...
            meanx(erpdata(:,:,pp.cond_diff{2},s),[1 2]);
        gfp_plot_data = var( meanx(erpdata(:,:,pp.cond_diff{1},s),[1 2]) - ...
            meanx(erpdata(:,:,pp.cond_diff{2},s),[1 2]) , 0, 2 );
        erp_topo_data = meanx(erpdata(t_start_ck:t_end_ck,:,pp.cond_diff{1},s),2) - ...
            meanx(erpdata(t_start_ck:t_end_ck,:,pp.cond_diff{2},s),2);
        ersp_plot_data = 10 * log10( bsxfun( @rdivide, squeeze(wave_totpowdata(:,ck_chan,:,pp.cond_diff{1},s)), ...
            squeeze(wave_totpowdata(:,ck_chan,:,pp.cond_diff{2},s)) ) );
        itc_plot_data = bsxfun( @minus, squeeze(wave_evknormdata(:,ck_chan,:,pp.cond_diff{1},s)), ...
            squeeze(wave_evknormdata(:,ck_chan,:,pp.cond_diff{2},s)) );
        ercoh_plot_data = bsxfun( @minus, meanx(cohdata(:,:,pp.cond_diff{1},plot_hypinds,s),[1 2]), ...
            meanx(cohdata(:,:,pp.cond_diff{2},plot_hypinds,s),[1 2]) );
    else
        erp_plot_data = squeeze(erpdata(:,:,cond,s));
        gfp_plot_data = var( meanx(erpdata(:,:,cond,s),[1 2]), 0, 2 );
        erp_topo_data = meanx(erpdata(t_start_ck:t_end_ck,:,cond,s),2);
        ersp_plot_data = 10 * log10( bsxfun( @rdivide, squeeze(wave_totpowdata(:,ck_chan,:,cond,s)), ...
            meanx(wave_totpowdata(t_start_b:t_end_b,ck_chan,:,:,s),3)' ) ); %dB-normed to common baseline
        itc_plot_data = bsxfun( @minus, squeeze(wave_evknormdata(:,ck_chan,:,cond,s)), ...
            meanx(wave_evknormdata(t_start_b:t_end_b,ck_chan,:,:,s),3)' ); %subtracting common baseline
        ercoh_plot_data = bsxfun( @minus, meanx(cohdata(:,:,cond,plot_hypinds,s),[1 2]), ...
            meanx(cohdata(t_start_b:t_end_b,:,:,plot_hypinds,s),2) ); %subtracting common baseline
    end
    %butterfly
    subplot_dummy=subplot_dummy+1; subplot(length(pp.plotn_cond),6,subplot_dummy);
    plot(erp_plot_data); %times are row, channels columns
    axis([scl.t_start_erp scl.t_end_erp -.8 .8]);
    vline(scl.t_zero,'k--'); hold off;
    set(gca,'XTick',scl.t_xtick_erp,'XTickLabel',scl.t_xtick_ms)

    % GFP / topographical variance
    subplot_dummy=subplot_dummy+1; subplot(length(pp.plotn_cond),6,subplot_dummy);
    plot(gfp_plot_data);
    axis([scl.t_start_erp scl.t_end_erp 0 .1]);
    vline(scl.t_zero,'k--'); hold off;
    set(gca,'XTick',scl.t_xtick_erp,'XTickLabel',scl.t_xtick_ms)
    
    %topography
    subplot_dummy=subplot_dummy+1; subplot(length(pp.plotn_cond),6,subplot_dummy);
    h=topoplot(erp_topo_data,chan_locs,'maplimits','absmax', ...
        'electrodes',pp.topo_elecs,'colormap',cmap,'style','fill','numcontour',pp.n_contour);
    set(h,'EdgeColor','None');
    
    %ersp
    subplot_dummy=subplot_dummy+1; subplot(length(pp.plotn_cond),6,subplot_dummy);
    [~,h]=contourf(fliplr(ersp_plot_data)',pp.n_contour);
    set(h,'EdgeColor','None'); axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms);
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label);
    grid on; set(gca,'Layer','Top'); hold on;
    plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
    colorbar;
    
    %itc
    subplot_dummy=subplot_dummy+1; subplot(length(pp.plotn_cond),6,subplot_dummy);
    [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
    set(h,'EdgeColor','None'); axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms);
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label);
    grid on; set(gca,'Layer','Top'); hold on;
    plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
    colorbar;
    
    %ercoh
    subplot_dummy=subplot_dummy+1; subplot(length(pp.plotn_cond),6,subplot_dummy);
    [~,h]=contourf(fliplr(ercoh_plot_data)',pp.n_contour);
    set(h,'EdgeColor','None'); axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms);
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label);
    grid on; set(gca,'Layer','Top'); hold on;
    plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
    colorbar;
    
end

adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,num2str(s),subplot_dims);

yn=input(sprintf('S%d:',s),'s');
if strcmpi(yn,'n')
    rejmat(s)=true;
else
    rejmat(s)=false;
end

    
end


%% plot topography of ERPs

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

erp_topo_scale=[-8 11];
erp_diff_limits=[-4 11];
cmap=makecmap(erp_topo_scale);
cmap_diff=makecmap(erp_diff_limits);

v=zeros(length(pp.plotn_g),length(pp.plotn_cond),length(pp.t_start_ms),2);
for group=pp.chosen_g(pp.plotn_g)
figure; subplot_dummy=0;
overtitle=sprintf('Topography of ERPs / %s',scl.g_label{group});
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)
    [~,t_start]=min(abs(scl.t_ms_erp-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms_erp-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    sp(cond)=subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy);
    %if win~=pp.maxwin
    if cond==imp.maxconds+1
        erp_topo_data = mean(mean(mean(erpdata(t_start:t_end,:,pp.cond_diff{1},s_inds_g(:,group)),1),3),4) - ...
    mean(mean(mean(erpdata(t_start:t_end,:,pp.cond_diff{2},s_inds_g(:,group)),1),3),4);
        h=topoplot(erp_topo_data,chan_locs,'maplimits',[erp_diff_limits(1) erp_diff_limits(2)], ...
            'electrodes',pp.topo_elecs,'colormap',cmap_diff,'style','fill','numcontour',pp.n_contour);
        %freezeColors;
    else
        erp_topo_data=mean(mean(erpdata(t_start:t_end,:,cond,s_inds_g(:,group)),1),4);
        h=topoplot(erp_topo_data,chan_locs,'maplimits',[erp_topo_scale(1) erp_topo_scale(2)], ...
            'electrodes',pp.topo_elecs,'colormap',cmap,'style','fill','numcontour',pp.n_contour);
        %freezeColors;
    end
    %shading flat
    %set(h,'EdgeColor','None');
    %else
    %    set(gca,'Visible','Off');
    %end
    %title(sprintf('%s, %d - %d ms, %s',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),scl.g_label{group}))
    v(group,cond,win,1)=min(erp_topo_data); v(group,cond,win,2)=max(erp_topo_data);
    end
    %color bar
    %{
    cb_pos = get(gca,'Position') + [0.04 0 0 0];
    cb_pos(3) = 0.05; %set width
    if cond==imp.maxconds+1
        cb_lim = erp_diff_limits;
    else
        cb_lim = erp_topo_scale;
    end
    colorscale([1 256], cb_lim, range(cb_lim)/5, 'vert', 2, ...
        'Position',cb_pos);
    freezeColors;
    %}
end
%adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims,true);
tightfig;
set_print_size(17,17/scl.phi);
end
c(1)=min(min(min(v(:,1:end-1,:,1)))); c(2)=max(max(max(v(:,1:end-1,:,2))));
c_diff(1)=min(min(v(:,end,:,1))); c_diff(2)=max(max(v(:,end,:,2)));

figure;
set(gcf,'position',[2800 120 300 400]);
colorscale_plot(erp_topo_scale, cmap, 0.25);
colorscale_plot(erp_diff_limits, cmap_diff, 0.75);

clear_plotassistvars

%% plot topography of ERPs from erpdata_h1

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

erp_topo_scale=[-8 11];
erp_diff_limits=[-4 11];
cmap=makecmap(erp_topo_scale);
cmap_diff=makecmap(erp_diff_limits);

v=zeros(length(pp.plotn_g),length(pp.plotn_cond),length(pp.t_start_ms),2);
for group=pp.chosen_g(pp.plotn_g)
figure; subplot_dummy=0;
overtitle=sprintf('Topography of ERPs / %s',scl.g_label{group});
for cond=pp_h1.plotn_cond
    for win=1:length(pp.t_start_ms)
    [~,t_start]=min(abs(scl_h1.t_ms_erp-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl_h1.t_ms_erp-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    sp(cond)=subplot(length(pp_h1.plotn_cond),length(pp.t_start_ms),subplot_dummy);
    if cond==imp.maxconds+1
        erp_topo_data = mean(mean(mean(erpdata(t_start:t_end,:,pp.cond_diff{1},s_inds_g(:,group)),1),3),4) - ...
    mean(mean(mean(erpdata(t_start:t_end,:,pp.cond_diff{2},s_inds_g(:,group)),1),3),4);
        h=topoplot(erp_topo_data,chan_locs,'maplimits',[erp_diff_limits(1) erp_diff_limits(2)], ...
            'electrodes',pp.topo_elecs,'colormap',cmap_diff,'style','fill','numcontour',pp.n_contour);
    else
        erp_topo_data=mean(mean(erpdata_h1(t_start:t_end,:,cond,s_inds_g(:,group)),1),4);
        h=topoplot(erp_topo_data,chan_locs,'maplimits',[erp_topo_scale(1) erp_topo_scale(2)], ...
            'electrodes',pp.topo_elecs,'colormap',cmap,'style','fill','numcontour',pp.n_contour);
    end
    v(group,cond,win,1)=min(erp_topo_data); v(group,cond,win,2)=max(erp_topo_data);
    end
end
tightfig;
set_print_size(17,17/scl.phi);
end
c(1)=min(min(min(v(:,1:end-1,:,1)))); c(2)=max(max(max(v(:,1:end-1,:,2))));
c_diff(1)=min(min(v(:,end,:,1))); c_diff(2)=max(max(v(:,end,:,2)));

figure;
set(gcf,'position',[2800 120 300 400]);
colorscale_plot(erp_topo_scale, cmap, 0.25);
colorscale_plot(erp_diff_limits, cmap_diff, 0.75);

clear_plotassistvars


%% plot topography of ERPs (HBNL TOPO VERSION)

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

erp_topo_scale=[-8 12];
erp_diff_limits=[-4 3];
v=zeros(length(pp.plotn_g),length(pp.plotn_cond),length(pp.t_start_ms),2);
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
        topoplot_hbnl(erp_topo_data','',erp_diff_limits(1),erp_diff_limits(2),...
            5, 2, 1, 0, pp.cmap, false);
    else
        erp_topo_data=mean(mean(erpdata(t_start:t_end,:,cond,s_inds_g(:,group)),1),4);
        topoplot_hbnl(erp_topo_data','',erp_topo_scale(1),erp_topo_scale(2),...
            15, 2, 1, 0, pp.cmap, false);
    end
    %title(sprintf('%s, %d - %d ms, %s',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),scl.g_label{group}))
    v(group,cond,win,1)=min(erp_topo_data); v(group,cond,win,2)=max(erp_topo_data);
    end
    colorbar;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
tightfig;
end
c(1)=min(min(min(v(:,1:end-1,:,1)))); c(2)=max(max(max(v(:,1:end-1,:,2))));
c_diff(1)=min(min(v(:,end,:,1))); c_diff(2)=max(max(v(:,end,:,2)));
clear_plotassistvars

%% look at EROs in individual bands / superimpose groups

ero_yax=[35 85];
erodiff_yax=[-4 14];

%h_line=zeros(length(pp.plotn_cond),length(pp.chosen_g));
for chan=pp.chosen_chan(pp.plotn_chan)
for freq_range=pp.plotn_f
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    figure;
    subplot_dummy=0;
    for cond=pp.plotn_cond
        subplot_dummy=subplot_dummy+1;
        sp(cond)=subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy);
        for group=pp.chosen_g(pp.plotn_g)
            if cond==imp.maxconds+1
                ero_plot_data = mean(mean(mean(wave_totpowdata(:,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),3),4),5) - ...
            mean(mean(mean(wave_totpowdata(:,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),3),4),5);
                ero_plot_data_std = std(mean(mean(wave_totpowdata(:,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),3),4) - ...
            mean(mean(wave_totpowdata(:,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),3),4),0,5)/sqrt(sum(s_inds_g(:,group)));
            else
                ero_plot_data=mean(mean(wave_totpowdata(:,chan,f_end:f_start,cond,s_inds_g(:,group)),3),5);
                ero_plot_data_std=std(mean(wave_totpowdata(:,chan,f_end:f_start,cond,s_inds_g(:,group)),3),0,5)/sqrt(sum(s_inds_g(:,group)));
            end
            %plot(ero_plot_data,scl.g_color{group}); hold on; %,'Color',scl.s_color(pp.chosen_s,:)); hold on
            h=shadedErrorBar(1:imp.maxtimepts,ero_plot_data,ero_plot_data_std, scl.g_color{group}); hold on;
            %h_line(cond,group)=h.mainLine;
        end
        if cond==imp.maxconds+1
            axis([1 scl.t_end erodiff_yax(1) erodiff_yax(2)]);
        else
            axis([1 scl.t_end ero_yax(1) ero_yax(2)]);
        end
        %axis tight
        vline(scl.t_zero,'k--'); hold off;
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
        grid on
        title([scl.chan_label{chan},'/',scl.cond_label{cond},'/',num2str(pp.f_start_hz(freq_range)),'-',num2str(pp.f_end_hz(freq_range)),' Hz'])
    end
tightfig;
linkaxes(sp(1:end-1))
end
end
clear_plotassistvars

%% look at EROS in individual bands / superimpose conditions

ero_yax=[-10 85];

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
                ero_plot_data = mean(mean(mean(wave_totpowdata(:,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),3),4),5) - ...
            mean(mean(mean(wave_totpowdata(:,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),3),4),5);
                ero_plot_data_std = std(mean(mean(wave_totpowdata(:,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),3),4) - ...
            mean(mean(wave_totpowdata(:,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),3),4),0,5)/sqrt(sum(s_inds_g(:,group)));
            else
                ero_plot_data=mean(mean(wave_totpowdata(:,chan,f_end:f_start,cond,s_inds_g(:,group)),3),5);
                ero_plot_data_std=std(mean(wave_totpowdata(:,chan,f_end:f_start,cond,s_inds_g(:,group)),3),0,5)/sqrt(sum(s_inds_g(:,group)));
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

sp_rowlabel=[];
%sp_columnlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';
subplot_dims=pp.sp_d;

pp.figdum=pp.figdum_init;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
for group=pp.chosen_g(pp.plotn_g)
for chan=25 %pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
clear overtitle
overtitle{pp.figdum}=sprintf('%s / %s',scl.chan_label{chan},scl.g_label{group});
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        ero_plot_data=squeeze(mean(mean(wave_totpowdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),4),5)-...
            mean(mean(wave_totpowdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),4),5));
    else
        ero_plot_data=squeeze(mean(wave_totpowdata(:,chan,:,cond,s_inds_g(:,group)),5));
    end
    [~,h]=contourf(fliplr(ero_plot_data)',pp.n_contour);
    %shading flat; 
    set(h,'EdgeColor','None');
    %imagesc(ero_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,cond,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
end
end
%
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
cmap_diff=makecmap(c_diff);
cmap=makecmap(c);
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[120 120 1500 400]);
for splot=1:subplot_dummy
    temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
    if splot==subplot_dummy
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
tightfig;
set_print_size(20,8);
plottitle(overtitle{fig});
end
%make color bars separately
figure;
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
clear_plotassistvars;

%% image ERO (ERSP) in time-freq at a chosen channel, in decibel change from baseline
% uses condition-mean (common) baseline

sp_rowlabel=[];
sp_columnlabel={scl.cond_label{pp.plotn_cond}};
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

pp.figdum=pp.figdum_init;
%pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for chan=[7 25] %pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
overtitle{pp.figdum}=sprintf('%s / %s',scl.chan_label{chan},scl.g_label{group});
ero_plot_base=squeeze( repmat( mean(mean(mean(wave_totpowdata(t_start_b:t_end_b,chan,:,:,s_inds_g(:,group)),1),4),5) , ...
    [1 length(scl.t_ms) 1] ) );
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        ero_plot_data = 10 * log10( ...
            squeeze(mean(wave_totpowdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),5)) ./ ...
            squeeze(mean(wave_totpowdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),5)) ); 
    else
        ero_plot_data = 10 * log10( ...
            squeeze(mean(wave_totpowdata(:,chan,:,cond,s_inds_g(:,group)),5)) ./ ... 
            ero_plot_base ); 
    end
    [~,h]=contourf(fliplr(ero_plot_data)',pp.n_contour);
    %shading flat; 
    set(h,'EdgeColor','None');
    %imagesc(ero_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,cond,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,pp.sp_d);
end
end
%
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
%c=symmetrize(c);
%c_diff=symmetrize(c_diff);
cmap_diff=makecmap(c_diff);
cmap=makecmap(c);
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[120 120 1500 400]);
for splot=1:subplot_dummy
    temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
    if splot==subplot_dummy
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
tightfig;
set_print_size(20,8);
plottitle(overtitle{fig},2);
end
%make color bars separately
figure;
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
clear_plotassistvars;

%% image ERO (ERSP) in time-freq at a chosen channel, in decibel change from baseline
% uses condition-mean (common) baseline
% AVERAGE AFTER dB calculation!

sp_rowlabel=[];
sp_columnlabel={scl.cond_label{pp.plotn_cond}};
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

p_alpha=.01;
%n_perms=round(1/p_alpha);
n_perms=500;
do_statmask=false;

pp.figdum=pp.figdum_init;
%pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for chan=7 %pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
overtitle{pp.figdum}=sprintf('%s / %s',scl.chan_label{chan},scl.g_label{group});
ero_plot_base=repmat(reshape(mean(mean(wave_totpowdata(t_start_b:t_end_b,chan,:,:,s_inds_g(:,group)),1),4), ...
    [1 length(scl.freqs) sum(s_inds_g(:,group))]),[length(scl.t_ms) 1 1]);
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        ero_plot_data = squeeze(mean( 10 * log10(  ...
            squeeze(wave_totpowdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group))) ./ ...
            squeeze(wave_totpowdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)))  ), 3) ); 
        if do_statmask
            [stats, df, pvals, surrog] = statcond( ero_diff_data, 'paired','on','method','perm', ...
                'naccu', n_perms, 'alpha', p_alpha, 'structoutput', 'on');
            ero_plot_data(~stats.mask) = 0;
        end
        imagesc(fliplr(ero_plot_data)');
        set(gca,'YDir','normal');
    else
        ero_plot_data = squeeze(mean( 10 * log10(  ...
            squeeze(wave_totpowdata(:,chan,:,cond,s_inds_g(:,group))) ./ ... 
            ero_plot_base), 3) ); 
        ero_diff_data{cond} = 10 * log10(meanx( ...
            wave_totpowdata(:,chan,:,cond,s_inds_g(:,group)),[1 3 5]));
        %ero_stat_data{1}= 10 * log10(bsxfun(@rdivide, ...
        %    squeeze(wave_totpowdata(:,chan,:,cond,s_inds_g(:,group))), ...
        %    mean(squeeze(wave_totpowdata(t_start_b:t_end_b,chan,:,cond,s_inds_g(:,group))),1) ...
        %    ));
        %ero_stat_data{2}=zeros(size(ero_stat_data{1}));
        %[stats, df, pvals, surrog] = statcond( ero_stat_data, 'paired','on','method','perm', ...
        %    'naccu', n_perms, 'alpha', p_alpha, 'structoutput', 'on');
        %ero_plot_data(~stats.mask) = 0;
        [~,h]=contourf(fliplr(ero_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
    end
    %shading flat; 
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,cond,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,pp.sp_d);
end
end
%
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
%c=symmetrize(c);
%c_diff=symmetrize(c_diff);
cmap=makecmap(c,0,'redblue',256);
cmap_diff=makecmap(c_diff,0,'redblue',256);
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[1350 120 1500 400]);
for splot=1:subplot_dummy
    temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
    if splot==subplot_dummy
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
%tightfig;
set_print_size(20,8);
plottitle(overtitle{fig},2);
end
%make color bars separately
figure;
set(gcf,'position',[2850 120 300 400]);
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
clear_plotassistvars;

%% image ERO POWER (ERSP) in time-freq at a chosen channel, in decibel change from baseline
% uses condition-mean (common) baseline
% AVERAGE AFTER dB calculation!

sp_rowlabel=[];
sp_columnlabel={scl.cond_label{pp.plotn_cond}};
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

p_alpha=.01;
%n_perms=round(1/p_alpha);
n_perms=500;
do_statmask=false;

pp.figdum=pp.figdum_init;
%pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for chan=[7] %pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
overtitle{pp.figdum}=sprintf('%s / %s',scl.chan_label{chan},scl.g_label{group});
ero_plot_base=repmat(reshape(mean(mean(wave_totpowdata(t_start_b:t_end_b,chan,:,:,s_inds_g(:,group)),1),4), ...
    [1 length(scl.freqs) sum(s_inds_g(:,group))]),[length(scl.t_ms) 1 1]);
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        ero_plot_data = squeeze(mean( 10 * log10(  ...
            squeeze(wave_totpowdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group))) ./ ...
            squeeze(wave_totpowdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)))  ), 3) ); 
        if do_statmask
            [stats, df, pvals, surrog] = statcond( ero_diff_data, 'paired','on','method','perm', ...
                'naccu', n_perms, 'alpha', p_alpha, 'structoutput', 'on');
            ero_plot_data(~stats.mask) = 0;
        end
        imagesc(fliplr(ero_plot_data)');
        set(gca,'YDir','normal');
    else
        ero_plot_data = squeeze(mean( 10 * log10(  ...
            squeeze(wave_totpowdata(:,chan,:,cond,s_inds_g(:,group))) ./ ... 
            ero_plot_base), 3) ); 
        ero_diff_data{cond} = 10 * log10(meanx( ...
            wave_totpowdata(:,chan,:,cond,s_inds_g(:,group)),[1 3 5]));
        %ero_stat_data{1}= 10 * log10(bsxfun(@rdivide, ...
        %    squeeze(wave_totpowdata(:,chan,:,cond,s_inds_g(:,group))), ...
        %    mean(squeeze(wave_totpowdata(t_start_b:t_end_b,chan,:,cond,s_inds_g(:,group))),1) ...
        %    ));
        %ero_stat_data{2}=zeros(size(ero_stat_data{1}));
        %[stats, df, pvals, surrog] = statcond( ero_stat_data, 'paired','on','method','perm', ...
        %    'naccu', n_perms, 'alpha', p_alpha, 'structoutput', 'on');
        %ero_plot_data(~stats.mask) = 0;
        [~,h]=contourf(fliplr(ero_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
    end
    %shading flat; 
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,cond,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,pp.sp_d);
end
end
%
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
%c=symmetrize(c);
%c_diff=symmetrize(c_diff);
cmap=makecmap(c,0,'redblue',256);
cmap_diff=makecmap(c_diff,0,'redblue',256);
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[1350 120 1500 400]);
for splot=1:subplot_dummy
    temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
    if splot==subplot_dummy
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
%tightfig;
set_print_size(20,8);
plottitle(overtitle{fig},2);
end
%make color bars separately
figure;
set(gcf,'position',[2850 120 300 400]);
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
clear_plotassistvars;

%% image ERO POWER (ERSP) in time-freq in a REGION, in decibel change from baseline
% uses condition-mean (common) baseline
% AVERAGE AFTER dB calculation!

sp_rowlabel=[];
sp_columnlabel={scl.cond_label{pp.plotn_cond}};
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

p_alpha=.01;
%n_perms=round(1/p_alpha);
n_perms=500;

%pp.figdum=pp.figdum_init;
pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for hyp=1:3
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
overtitle{pp.figdum}=sprintf('%s / %s',opt.pair_indlbls{hyp},scl.g_label{group});
plot_hypinds=find(opt.pair_inds==hyp);
plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
ero_plot_base=mean( mean( wave_totpowdata(t_start_b:t_end_b,plot_hypchans,:,:,s_inds_g(:,group)), 1), 4);
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        ero_plot_data = 10*log10( meanx( bsxfun(@rdivide, ...
            wave_totpowdata(:,plot_hypchans,:,pp.cond_diff{1},s_inds_g(:,group)), ...
            wave_totpowdata(:,plot_hypchans,:,pp.cond_diff{2},s_inds_g(:,group))), [1 3]));
        %[stats, df, pvals, surrog] = statcond( ero_diff_data, 'paired','on','method','perm', ...
        %    'naccu', n_perms, 'alpha', p_alpha, 'structoutput', 'on');
        %ero_plot_data(~stats.mask) = 0;
        imagesc(fliplr(ero_plot_data)');
        set(gca,'YDir','normal');
    else
        ero_plot_data = 10*log10( meanx( bsxfun(@rdivide, ...
            wave_totpowdata(:,plot_hypchans,:,cond,s_inds_g(:,group)), ...
            ero_plot_base), [1 3]));
        ero_diff_data{cond} = 10 * log10(meanx( ...
            wave_totpowdata(:,plot_hypchans,:,cond,s_inds_g(:,group)),[1 3 5]));
        [~,h]=contourf(fliplr(ero_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
    end
    %shading flat; 
    %imagesc(ero_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,hyp,cond,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,pp.sp_d);
end
end
%
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
cmap=makecmap(c);
cmap_diff=makecmap(c_diff,0,'redblue');
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[120 120 1500 400]);
for splot=1:subplot_dummy
    temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
    if splot==subplot_dummy
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
tightfig;
set_print_size(20,8);
plottitle(overtitle{fig});
end
%make color bars separately
figure;
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
clear_plotassistvars;


%% image ERO (ERSP) in time-freq in a REGION, in decibel change from baseline
% uses condition-mean (common) baseline
% AVERAGE AFTER dB calculation!

sp_rowlabel=[];
sp_columnlabel={scl.cond_label{pp.plotn_cond}};
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

p_alpha=.01;
%n_perms=round(1/p_alpha);
n_perms=500;

%pp.figdum=pp.figdum_init;
pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for hyp=1:3
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
overtitle{pp.figdum}=sprintf('%s / %s',opt.pair_indlbls{hyp},scl.g_label{group});
plot_hypinds=find(opt.pair_inds==hyp);
plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
ero_plot_base=mean( mean( wave_totpowdata(t_start_b:t_end_b,plot_hypchans,:,:,s_inds_g(:,group)), 1), 4);
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        ero_plot_data = 10*log10( meanx( bsxfun(@rdivide, ...
            wave_totpowdata(:,plot_hypchans,:,pp.cond_diff{1},s_inds_g(:,group)), ...
            wave_totpowdata(:,plot_hypchans,:,pp.cond_diff{2},s_inds_g(:,group))), [1 3]));
        [stats, df, pvals, surrog] = statcond( ero_diff_data, 'paired','on','method','perm', ...
            'naccu', n_perms, 'alpha', p_alpha, 'structoutput', 'on');
        ero_plot_data(~stats.mask) = 0;
        imagesc(fliplr(ero_plot_data)');
        set(gca,'YDir','normal');
    else
        ero_plot_data = 10*log10( meanx( bsxfun(@rdivide, ...
            wave_totpowdata(:,plot_hypchans,:,cond,s_inds_g(:,group)), ...
            ero_plot_base), [1 3]));
        ero_diff_data{cond} = 10 * log10(meanx( ...
            wave_totpowdata(:,plot_hypchans,:,cond,s_inds_g(:,group)),[1 3 5]));
        [~,h]=contourf(fliplr(ero_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
    end
    %shading flat; 
    %imagesc(ero_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,hyp,cond,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,pp.sp_d);
end
end
%
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
cmap=makecmap(c);
cmap_diff=makecmap(c_diff,0,'redblue');
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[120 120 1500 400]);
for splot=1:subplot_dummy
    temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
    if splot==subplot_dummy
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
tightfig;
set_print_size(20,8);
plottitle(overtitle{fig});
end
%make color bars separately
figure;
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
clear_plotassistvars;

%% plot ERO as a topographic plot over the head

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

ero_topo_scale=[-.24 .21];
ero_diff_limits=[-0.02 .09];
cmap=makecmap(ero_topo_scale);
cmap_diff=makecmap(ero_diff_limits);

v=zeros(length(pp.plotn_cond),length(pp.chosen_g),length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms),2);
for group=pp.chosen_g(pp.plotn_g)
for freq_range=1:3 %pp.plotn_f
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
        sp_temp=subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy);
        if win~=pp.maxwin
        if cond==imp.maxconds+1
            topo_data=squeeze(mean(mean(mean(mean(wave_totpowdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4),5)) - ...
                squeeze(mean(mean(mean(mean(wave_totpowdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4),5));
            h=topoplot(topo_data,chan_locs,'maplimits',[ero_diff_limits(1) ero_diff_limits(2)], ...
                'electrodes',pp.topo_elecs,'colormap',cmap_diff,'style','fill','numcontour',7);
            colormap(sp_temp,cmap_diff);
        else
            %topo_data=squeeze(mean(mean(mean(wave_totpowdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5));
            topo_data=squeeze(mean(mean(mean(wave_totpowdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5)) -...
                squeeze(mean(mean(mean(wave_totpowdata(1:scl.t_zero,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5));
            h=topoplot(topo_data,chan_locs,'maplimits',[ero_topo_scale(1) ero_topo_scale(2)], ...
                'electrodes',pp.topo_elecs,'colormap',cmap,'style','fill','numcontour',7);
            colormap(sp_temp,cmap);
        end
        set(h,'EdgeColor','None');
        else
            set(gca,'Visible','Off');
        end
        v(cond,group,freq_range,win,1) = min(topo_data); v(cond,group,freq_range,win,2) = max(topo_data);        
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
tightfig;
set_print_size(17,17/scl.phi);
end
end
figure;
colorscale_plot(ero_topo_scale, cmap, 0.25);
colorscale_plot(ero_diff_limits, cmap_diff, 0.75);
c(1)=min(min(min(min(v(1:end-1,:,:,:,1))))); c(2)=max(max(max(max(v(1:end-1,:,:,:,2)))));
c_diff(1)=min(min(min(v(end,:,:,:,1)))); c_diff(2)=max(max(max(v(end,:,:,:,2))));
clear_plotassistvars

%% plot ERO as a topographic plot over the head (dB-normed to common baseline)

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

ero_topo_scale=[-2 4.2];
ero_diff_limits=[-0.8 .9];
cmap=makecmap(ero_topo_scale);
cmap_diff=makecmap(ero_diff_limits);

v=zeros(length(pp.plotn_cond),length(pp.chosen_g),length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms),2);
clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for freq_range=1:3
figure; subplot_dummy=0;
overtitle=sprintf('Topography of EROs in %1.1f - %1.1f Hz / %s',pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.g_label{group});
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
ero_topo_base=meanx(wave_totpowdata(t_start_b:t_end_b,:,f_end:f_start,:,s_inds_g(:,group)), 2);
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        subplot_dummy=subplot_dummy+1;
        sp_temp=subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy);
        if cond==imp.maxconds+1
            topo_data=10*log10( meanx(wave_totpowdata(t_start:t_end,:,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),2) ./ ...
                ero_topo_base ) - ...
                10*log10( meanx(wave_totpowdata(t_start:t_end,:,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),2) ./ ...
                ero_topo_base );
            h=topoplot(topo_data,chan_locs,'maplimits',[ero_diff_limits(1) ero_diff_limits(2)], ...
                'electrodes',pp.topo_elecs,'colormap',cmap_diff,'style','fill','numcontour',7);
            colormap(sp_temp,cmap_diff);
        else
            topo_data=10*log10( meanx(wave_totpowdata(t_start:t_end,:,f_end:f_start,cond,s_inds_g(:,group)),2) ./ ...
                ero_topo_base );
            h=topoplot(topo_data,chan_locs,'maplimits',[ero_topo_scale(1) ero_topo_scale(2)], ...
                'electrodes',pp.topo_elecs,'colormap',cmap,'style','fill','numcontour',7);
            colormap(sp_temp,cmap);
        end
        set(h,'EdgeColor','None');
        v(cond,group,freq_range,win,1) = min(topo_data); v(cond,group,freq_range,win,2) = max(topo_data);        
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
tightfig;
set_print_size(17,17/scl.phi);
end
end
c(1)=min(min(min(min(v(1:end-1,:,:,:,1))))); c(2)=max(max(max(max(v(1:end-1,:,:,:,2)))));
c_diff(1)=min(min(min(v(end,:,:,:,1)))); c_diff(2)=max(max(max(v(end,:,:,:,2))));
figure;
colorscale_plot(ero_topo_scale, cmap, 0.25);
colorscale_plot(ero_diff_limits, cmap_diff, 0.75);
clear_plotassistvars

%% plot ERO as a topographic plot over the head (dB-normed to common baseline)
% average after dB transformation

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

ero_topo_scale=[-2 4.2];
ero_diff_limits=[-0.8 .9];
cmap=makecmap(ero_topo_scale,0,'redblue');
cmap_diff=makecmap(ero_diff_limits,0,'redblue');

v=zeros(length(pp.plotn_cond),length(pp.chosen_g),length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms),2);
clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for freq_range=2
figure; subplot_dummy=0;
overtitle=sprintf('Topography of EROs in %1.1f - %1.1f Hz / %s',pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.g_label{group});
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
ero_topo_base=meanx(wave_totpowdata(t_start_b:t_end_b,:,f_end:f_start,:,s_inds_g(:,group)), [2 5]);
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        subplot_dummy=subplot_dummy+1;
        sp_temp=subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy);
        if win~=pp.maxwin
        if cond==imp.maxconds+1
            topo_data=squeeze(mean( 10*log10(meanx(wave_totpowdata(t_start:t_end,:,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),[2 5]) ./ ...
                meanx(wave_totpowdata(t_start:t_end,:,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),[2 5])) ,2));
            h=topoplot(topo_data,chan_locs,'maplimits',[ero_diff_limits(1) ero_diff_limits(2)], ...
                'electrodes',pp.topo_elecs,'colormap',cmap_diff,'style','fill','numcontour',7);
            colormap(sp_temp,cmap_diff);
        else
            topo_data=squeeze(mean( 10*log10(meanx(wave_totpowdata(t_start:t_end,:,f_end:f_start,cond,s_inds_g(:,group)),[2 5]) ./ ...
                ero_topo_base), 2));
            h=topoplot(topo_data,chan_locs,'maplimits',[ero_topo_scale(1) ero_topo_scale(2)], ...
                'electrodes',pp.topo_elecs,'colormap',cmap,'style','fill','numcontour',7);
            colormap(sp_temp,cmap);
        end
        set(h,'EdgeColor','None');
        else
            set(gca,'Visible','Off');
        end
        v(cond,group,freq_range,win,1) = min(topo_data); v(cond,group,freq_range,win,2) = max(topo_data);        
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
tightfig;
set_print_size(17,17/scl.phi);
end
end
c(1)=min(min(min(min(v(1:end-1,:,:,:,1))))); c(2)=max(max(max(max(v(1:end-1,:,:,:,2)))));
c_diff(1)=min(min(min(v(end,:,:,:,1)))); c_diff(2)=max(max(max(v(end,:,:,:,2))));
figure;
colorscale_plot(ero_topo_scale, cmap, 0.25);
colorscale_plot(ero_diff_limits, cmap_diff, 0.75);
clear_plotassistvars

%% plot ERO power as a topographic plot over the head (dB-normed to common baseline)
% average after dB transformation

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

ero_topo_scale=[-2 4.2];
ero_diff_limits=[-0.8 .9];
cmap=makecmap(ero_topo_scale,0,'redblue');
cmap_diff=makecmap(ero_diff_limits,0,'redblue');

v=zeros(length(pp.plotn_cond),length(pp.chosen_g),length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms),2);
clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for freq_range=2
figure; subplot_dummy=0;
overtitle=sprintf('Topography of EROs in %1.1f - %1.1f Hz / %s',pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.g_label{group});
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
ero_topo_base=meanx(wave_totpowdata(t_start_b:t_end_b,:,f_end:f_start,:,s_inds_g(:,group)), [2 5]);
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        subplot_dummy=subplot_dummy+1;
        sp_temp=subplot(length(pp.plotn_cond),length(pp.t_start_ms),subplot_dummy);
        if win~=pp.maxwin
        if cond==imp.maxconds+1
            topo_data=squeeze(mean( 10*log10(meanx(wave_totpowdata(t_start:t_end,:,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),[2 5]) ./ ...
                meanx(wave_totpowdata(t_start:t_end,:,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),[2 5])) ,2));
            h=topoplot(topo_data,chan_locs,'maplimits',[ero_diff_limits(1) ero_diff_limits(2)], ...
                'electrodes',pp.topo_elecs,'colormap',cmap_diff,'style','fill','numcontour',7);
            colormap(sp_temp,cmap_diff);
        else
            topo_data=squeeze(mean( 10*log10(meanx(wave_totpowdata(t_start:t_end,:,f_end:f_start,cond,s_inds_g(:,group)),[2 5]) ./ ...
                ero_topo_base), 2));
            h=topoplot(topo_data,chan_locs,'maplimits',[ero_topo_scale(1) ero_topo_scale(2)], ...
                'electrodes',pp.topo_elecs,'colormap',cmap,'style','fill','numcontour',7);
            colormap(sp_temp,cmap);
        end
        set(h,'EdgeColor','None');
        else
            set(gca,'Visible','Off');
        end
        v(cond,group,freq_range,win,1) = min(topo_data); v(cond,group,freq_range,win,2) = max(topo_data);        
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
tightfig;
set_print_size(17,17/scl.phi);
end
end
c(1)=min(min(min(min(v(1:end-1,:,:,:,1))))); c(2)=max(max(max(max(v(1:end-1,:,:,:,2)))));
c_diff(1)=min(min(min(v(end,:,:,:,1)))); c_diff(2)=max(max(max(v(end,:,:,:,2))));
figure;
colorscale_plot(ero_topo_scale, cmap, 0.25);
colorscale_plot(ero_diff_limits, cmap_diff, 0.75);
clear_plotassistvars


%% plot ERO as a topographic plot over the head (HBNL TOPOPLOT VERSION)

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

ero_topo_scale=[-2 18];
ero_diff_limits=[-9 3];

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
            topo_data=squeeze(mean(mean(mean(mean(wave_totpowdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4),5)) - ...
                squeeze(mean(mean(mean(mean(wave_totpowdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4),5));
            topoplot(topo_data,chan_locs,'maplimits',[ero_diff_limits(1) ero_diff_limits(2)],'electrodes',pp.topo_elecs,'colormap',pp.cmap);
        else
            %topo_data=squeeze(mean(mean(mean(wave_totpowdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5));
            topo_data=squeeze(mean(mean(mean(wave_totpowdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5)) -...
                squeeze(mean(mean(mean(wave_totpowdata(1:scl.t_zero,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5));
            topoplot(topo_data,chan_locs,'maplimits',[ero_topo_scale(1) ero_topo_scale(2)],'electrodes',pp.topo_elecs,'colormap',pp.cmap);
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

%% scatter ERO with behavioral measures

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
y_plotlabel='ERSP';

%for po=1:4

%x_plotlabel='Average Bet';
%x_plotlabel='Criterion';
%x_plotlabel=['Criterion Following ',behdata.cond1{po}];
x_plotlabel='Age';
%behav_data=behdata.avgbet;
%behav_data=behdata.crit_po(po,:);
behav_data=s_demogs.age_eeg';
%behav_data=behdata.crit;
%behav_data=behdata.avgbet_po(po,:); %avgbet; %

ero_axes=[min(min(behav_data)) max(max(behav_data)) -2.4 3.5]; %max(max(behav_data))
ero_axes_diff=[min(min(behav_data)) max(max(behav_data)) -1.2 2];

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_start_hz(pp.plotn_f)),pp.maxwin);
v=zeros(length(pp.plotn_chan),length(pp.plotn_f),length(pp.plotn_cond),pp.maxwin,2);
for chan=[7] %pp.chosen_chan(pp.plotn_chan)
for freq_range=1:2 %pp.plotn_f(6)
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    figure; subplot_dummy=0;
    overtitle=sprintf('ERSP at %s in %1.1f - %1.1f Hz vs. Behavior', ...
        scl.chan_label{chan},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
    for cond=pp.plotn_cond
    for win=1:pp.maxwin
    [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
    %    
    for group=pp.chosen_g(pp.plotn_g)
        
    ero_scatterdata_x=behav_data(s_inds_g(:,group))';
    ero_scatterbase = permute(meanx( ...
        wave_totpowdata(t_start_b:t_end_b,chan,f_end:f_start,:,s_inds_g(:,group)), ...
        [3 5]),[3 1 2]);
    
    if cond==imp.maxconds+1
        ero_scatterdata_y= meanx( 10*log10( bsxfun(@rdivide, ...
            squeeze(wave_totpowdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group))), ...
            squeeze(wave_totpowdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group))) ) ) ,3); 
        %[p_tfwin(:,chan,freq_range,win) stats]=robustfit(ero_scatterdata_x,ero_scatterdata_y);
        [p_tfwin(:,chan,freq_range,win),~,~,~,stats]=regress(ero_scatterdata_y,[ones(length(ero_scatterdata_x),1) ero_scatterdata_x]);
        scatter_h(group)=scatter(ero_scatterdata_x,ero_scatterdata_y, scl.g_color{group}); hold on;
        %plot(linspace(ero_axes_diff(1),ero_axes_diff(2),100),linspace(ero_axes_diff(3),ero_axes_diff(4),100),'k--'); hold on;
        plot(linspace(ero_axes_diff(1),ero_axes_diff(2),100),linspace(ero_axes_diff(1),ero_axes_diff(2),100)*...
            p_tfwin(2,chan,freq_range,win)+p_tfwin(1,chan,freq_range,win), scl.g_color{group}); hold on;
        if stats(3) < .05 %stats.p(2) < .05
            text((ero_axes_diff(1)+ero_axes_diff(2))/2,(ero_axes_diff(4))*(group/12)+1.5,['p=',num2str(stats(3),3)],'Color','k');
            text((ero_axes_diff(1)+ero_axes_diff(2))/2,(ero_axes_diff(4))*(group/12)+1,['r^2=',num2str(stats(1),3)],'Color','k');
        end
        axis(ero_axes_diff);
    else
        ero_scatterdata_y= meanx( 10*log10( bsxfun(@rdivide, ...
            squeeze(wave_totpowdata(t_start:t_end,chan,f_end:f_start,cond,s_inds_g(:,group))), ...
            ero_scatterbase ) ) ,3); 
        %[p_tfwin(:,chan,freq_range,win) stats]=robustfit(ero_scatterdata_x,ero_scatterdata_y);
        [p_tfwin(:,chan,freq_range,win),~,~,~,stats]=regress(ero_scatterdata_y,[ones(size(ero_scatterdata_x)) ero_scatterdata_x]);
        scatter_h(group)=scatter(ero_scatterdata_x,ero_scatterdata_y, scl.g_color{group}); hold on;
        %plot(linspace(ero_axes(1),ero_axes(2),100),linspace(ero_axes(3),ero_axes(4),100),'k--'); hold on;
        plot(linspace(ero_axes(1),ero_axes(2),100),linspace(ero_axes(1),ero_axes(2),100)*...
            p_tfwin(2,chan,freq_range,win)+p_tfwin(1,chan,freq_range,win), scl.g_color{group}); hold on;
        if stats(3) < .05 %stats.p(2) < .05
            text((ero_axes(1)+ero_axes(2))/2,(ero_axes(4))*(group/12)+2.5,['p=',num2str(stats(3),3)],'Color','k');
            text((ero_axes(1)+ero_axes(2))/2,(ero_axes(4))*(group/12)+2,['r^2=',num2str(stats(1),3)],'Color','k');
        end
        axis(ero_axes);
    end
    v(chan,freq_range,cond,win,1) = min(ero_scatterdata_y); v(chan,freq_range,cond,win,2) = max(ero_scatterdata_y);
    end
    hold off; grid on;
    end
    end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
tightfig; %dragzoom;
end
end
c(1)=min(min(min(min(v(:,:,1:end-1,:,1))))); c(2)=max(max(max(max(v(:,:,1:end-1,:,2)))));
c_diff(1)=min(min(min(v(:,:,end,:,1)))); c_diff(2)=max(max(max(v(:,:,end,:,2))));

%end

clear_plotassistvars

%% scatter ITC with behavioral measures

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
y_plotlabel='ITC';

%for po=1:4

%x_plotlabel='Average Bet';
%x_plotlabel=['Average Bet Following ',behdata.cond1{po}];
%behav_data=behdata.avgbet;
%behav_data=behdata.avgbet_po(po,:); %avgbet; %
%x_plotlabel=['Criterion Following ',behdata.cond1{po}];
%behav_data=behdata.crit_po(po,:);
x_plotlabel='Age';
behav_data=s_demogs.age_eeg';

itc_axes=[min(min(behav_data)) max(max(behav_data)) -.1 .4]; %max(max(behav_data))
itc_axes_diff=[min(min(behav_data)) max(max(behav_data)) -0.4 0.4];

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_start_hz(pp.plotn_f)),pp.maxwin);
v=zeros(length(pp.plotn_chan),length(pp.plotn_f),length(pp.plotn_cond),pp.maxwin,2);
for chan=[7] %pp.chosen_chan(pp.plotn_chan)
for freq_range=1:2 %pp.plotn_f(6)
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    figure; subplot_dummy=0;
    overtitle=sprintf('ITC at %s in %1.1f - %1.1f Hz vs. Behavior', ...
        scl.chan_label{chan},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
    for cond=pp.plotn_cond
    for win=1:pp.maxwin
    [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
    [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.plotn_cond),pp.maxwin,subplot_dummy);
    %    
    for group=pp.chosen_g(pp.plotn_g)
        
    itc_scatterdata_x=behav_data(s_inds_g(:,group))';
    itc_scatterbase = permute(meanx( ...
        wave_evknormdata(t_start_b:t_end_b,chan,f_end:f_start,:,s_inds_g(:,group)), ...
        [3 5]),[3 1 2]);
    
    if cond==imp.maxconds+1
        itc_scatterdata_y= meanx( bsxfun(@minus, ...
            squeeze(wave_evknormdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group))), ...
            squeeze(wave_evknormdata(t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group))) ) ,3); 
        %[p_tfwin(:,chan,freq_range,win) stats]=robustfit(itc_scatterdata_x,itc_scatterdata_y);
        [p_tfwin(:,chan,freq_range,win),~,~,~,stats]=regress(itc_scatterdata_y,[ones(length(itc_scatterdata_x),1) itc_scatterdata_x]);
        scatter_h(group)=scatter(itc_scatterdata_x,itc_scatterdata_y, scl.g_color{group}); hold on;
        %plot(linspace(itc_axes_diff(1),itc_axes_diff(2),100),linspace(itc_axes_diff(3),itc_axes_diff(4),100),'k--'); hold on;
        plot(linspace(itc_axes_diff(1),itc_axes_diff(2),100),linspace(itc_axes_diff(1),itc_axes_diff(2),100)*...
            p_tfwin(2,chan,freq_range,win)+p_tfwin(1,chan,freq_range,win), scl.g_color{group}); hold on;
        if stats(3) < .05
            text((itc_axes_diff(1)+itc_axes_diff(2))/2,(itc_axes_diff(4))*(group/12)+0.4,['p=',num2str(stats(3),3)],'Color','k');
            text((itc_axes_diff(1)+itc_axes_diff(2))/2,(itc_axes_diff(4))*(group/12)+0.3,['r^2=',num2str(stats(1),3)],'Color','k');
        end
        axis(itc_axes_diff);
    else
        itc_scatterdata_y= meanx( bsxfun(@minus, ...
            squeeze(wave_evknormdata(t_start:t_end,chan,f_end:f_start,cond,s_inds_g(:,group))), ...
            itc_scatterbase ) ,3); 
        %[p_tfwin(:,chan,freq_range,win) stats]=robustfit(itc_scatterdata_x,itc_scatterdata_y);
        [p_tfwin(:,chan,freq_range,win),~,~,~,stats]=regress(itc_scatterdata_y,[ones(length(itc_scatterdata_x),1) itc_scatterdata_x]);
        scatter_h(group)=scatter(itc_scatterdata_x,itc_scatterdata_y, scl.g_color{group}); hold on;
        %plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100),'k--'); hold on;
        plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(1),itc_axes(2),100)*...
            p_tfwin(2,chan,freq_range,win)+p_tfwin(1,chan,freq_range,win), scl.g_color{group}); hold on;
        if stats(3) < .05
            text((itc_axes(1)+itc_axes(2))/2,(itc_axes(4))*(group/12)+0.3,['p=',num2str(stats(3),3)],'Color','k');
            text((itc_axes(1)+itc_axes(2))/2,(itc_axes(4))*(group/12)+0.2,['r^2=',num2str(stats(1),3)],'Color','k');
        end
        axis(itc_axes);
    end
    v(chan,freq_range,cond,win,1) = min(itc_scatterdata_y); v(chan,freq_range,cond,win,2) = max(itc_scatterdata_y);
    end
    hold off; grid on;
    end
    end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle);
tightfig; %dragzoom;
end
end
c(1)=min(min(min(min(v(:,:,1:end-1,:,1))))); c(2)=max(max(max(max(v(:,:,1:end-1,:,2)))));
c_diff(1)=min(min(min(v(:,:,end,:,1)))); c_diff(2)=max(max(max(v(:,:,end,:,2))));

%end

clear_plotassistvars


%% image phase in time-freq at a chosen channel

v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
for group=pp.chosen_g
for chan=pp.chosen_chan(pp.plotn_chan)
figure; subplot_dummy=0;
for cond=1:imp.maxconds
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
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
    subplot(pp.sp_d(1),pp.sp_d(2),splot);
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
            topoplot(contour_wrapedges(topo_data),chan_locs,'maplimits',[phase_diff_limits(1) phase_diff_limits(2)],'electrodes',pp.topo_elecs,'colormap',pp.cmap);
        else
            topo_data=squeeze(angle(mean(mean(wave_evkdata(t_end:t_end,:,f_indiv,cond,s_inds_g(:,group)),1),5)));
            topoplot(contour_wrapedges(topo_data),chan_locs,'maplimits',[phase_topo_scale(1) phase_topo_scale(2)],'electrodes',pp.topo_elecs,'colormap',pp.cmap);
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

sp_rowlabel={''};
sp_columnlabel=scl.cond_label;
%sp_rowlabel=[];
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';
subplot_dims=pp.sp_d;

%pp.figdum=pp.figdum_init;
pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
for group=pp.chosen_g(pp.plotn_g)
for chan=[7] % 16 25 57 58] %pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
overtitle{pp.figdum}=sprintf('%s / %s',scl.chan_label{chan},scl.g_label{group});
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        itc_plot_data=squeeze(mean(mean(itcdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),4),5)-...
            mean(mean(itcdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),4),5));
        %itc_plot_data=squeeze(mean(wavelet_tot(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),5)-...
        %    mean(wavelet_tot(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),5));
    else
        itc_plot_data=squeeze(mean(itcdata(:,chan,:,cond,s_inds_g(:,group)),5));
        %itc_plot_data=squeeze(mean(wavelet_tot(:,chan,:,cond,s_inds_g(:,group)),5));
    end
    [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
    set(h,'EdgeColor','None');
    %imagesc(itc_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,subplot_dummy,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on;
    %title(['ITC at ',scl.chan_label{chan},' : ',scl.g_label{group},' : ',scl.cond_label{cond}])
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
end
end
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
cmap=makecmap(c);
cmap_diff=makecmap(c_diff);
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[120 120 1500 600]);
for splot=1:length(pp.plotn_cond);
    temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
    if splot==length(pp.plotn_cond)
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
tightfig;
set_print_size(20,8);
plottitle(overtitle{fig});
end
%make color bars separately
cb_pos=[0.25 0.1 0.05 0.8];
figure;
h=colorscale([1 256], c_diff, range(c_diff)/5, 'vert','Position',cb_pos);
colormap(h,cmap_diff);
cb_pos=[0.75 0.1 0.05 0.8];
h=colorscale([1 256], c, range(c)/5, 'vert','Position',cb_pos);
colormap(h,cmap);
%
clear_plotassistvars

%% image PURE ITC in time-freq at a chosen channel

sp_rowlabel={};
sp_columnlabel={};
%sp_rowlabel=[];
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';
subplot_dims=pp.sp_d;

%pp.figdum=pp.figdum_init;
pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
for group=pp.chosen_g(pp.plotn_g)
for chan=[25] % 16 25 57 58] %pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
overtitle{pp.figdum}=sprintf('%s / %s',scl.chan_label{chan},scl.g_label{group});
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        itc_plot_data=squeeze(mean(mean(wave_evknormdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),4),5)-...
            mean(mean(wave_evknormdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),4),5));
    else
        itc_plot_data=squeeze(mean(wave_evknormdata(:,chan,:,cond,s_inds_g(:,group)),5));
    end
    [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
    set(h,'EdgeColor','None');
    %shading flat;
    %imagesc(itc_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    v(group,chan,subplot_dummy,:) = caxis;
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    %title(['ITC at ',scl.chan_label{chan},' : ',scl.g_label{group},' : ',scl.cond_label{cond}])
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
end
end
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
cmap=makecmap(c);
cmap_diff=makecmap(c_diff);
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[120 120 1500 600]);
for splot=1:length(pp.plotn_cond);
    temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
    if splot==length(pp.plotn_cond)
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
tightfig;
set_print_size(20,8);
plottitle(overtitle{fig},2);
end
%make color bars separately
figure;
set(gcf,'position',[2800 120 300 400]);
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
clear_plotassistvars

%% image PURE ITC in time-freq at a chosen channel, as increase from baseline

sp_rowlabel={};
%sp_columnlabel=scl.cond_label;
sp_columnlabel={};
%sp_rowlabel=[];
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';
subplot_dims=pp.sp_d;

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

p_alpha=0.01;
n_perms=1000;
do_statmask=false;

pp.figdum=pp.figdum_init;
%pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
clear overtitle
for group=pp.chosen_g(pp.plotn_g)
for chan=7 %pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
overtitle{pp.figdum}=sprintf('%s / %s',scl.chan_label{chan},scl.g_label{group});
itc_plot_base=permute(meanx(wave_evknormdata(t_start_b:t_end_b,chan,:,:,s_inds_g(:,group)),[3 5]),[3 1 2]);
%itc_plot_base=zeros(size(permute(meanx(wave_evknormdata(t_start_b:t_end_b,chan,:,:,s_inds_g(:,group)),[3 5]),[3 1 2])));
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        itc_plot_data=squeeze(mean(mean(wave_evknormdata(:,chan,:,pp.cond_diff{1},s_inds_g(:,group)),4),5)-...
            mean(mean(wave_evknormdata(:,chan,:,pp.cond_diff{2},s_inds_g(:,group)),4),5));
        v(group,chan,cond,:) = [minx(itc_plot_data) maxx(itc_plot_data)];
        if do_statmask
            [stats, df, pvals, surrog] = statcond( itc_diff_data, 'paired','on','method','perm', ...
                'naccu', n_perms, 'alpha', p_alpha, 'structoutput', 'on');
            %itc_plot_data(~(stats.mask)) = 0;
            [p_fdr, p_masked] = fdr( stats.pval, p_alpha, 'parametric');
            %[p_masked, p_fdr, adj_ci_cvrg, adj_p] = fdr_bh(pvals, p_alpha, 'pdep');
            itc_plot_data(~(p_masked )) = 0;
        end
        imagesc(fliplr(itc_plot_data)');
        set(gca, 'YDir', 'normal');
        %itc_plot_conddiff_data = sign(itc_plot_data) .* -log10(stats.pval);
        %[~,h]=contourf(fliplr(itc_plot_conddiff_data)',pp.n_contour);
        %set(h,'EdgeColor','None');
        
    else
        %itc_diff_data{cond}=bsxfun(@minus,squeeze(wave_evknormdata(:,chan,:,cond,s_inds_g(:,group))), ...
        %    itc_plot_base);
        itc_diff_data{cond}=squeeze(wave_evknormdata(:,chan,:,cond,s_inds_g(:,group)));
        itc_plot_data = squeeze(mean( itc_diff_data{cond} , 3));
        [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
        v(group,chan,cond,:) = caxis;
    end
    
    %shading flat;
    %imagesc(itc_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    %title(['ITC at ',scl.chan_label{chan},' : ',scl.g_label{group},' : ',scl.cond_label{cond}])
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle{pp.figdum},subplot_dims);
end
end
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
cmap=makecmap(c, Inf, 'redgreen');
cmap_diff=makecmap(c_diff, 0, 'redgreen');
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[1300 120 1500 400]);
for splot=1:length(pp.plotn_cond);
    temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
    if splot==length(pp.plotn_cond)
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
%tightfig;
set_print_size(18,6);
plottitle(overtitle{fig},2);
end
%make color bars separately
cb_pos=[0.25 0.1 0.05 0.8];
figure;
set(gcf,'position',[2800 120 300 400]);
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
clear_plotassistvars

%% image ITC in time-freq in a REGION, as increase from baseline
% uses condition-mean (common) baseline

sp_rowlabel=[];
sp_columnlabel={scl.cond_label{pp.plotn_cond}};
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

p_alpha=0.01;
n_perms=2000;

%pp.figdum=pp.figdum_init;
pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),6,length(pp.plotn_cond),2);
clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for hyp=1:3
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
if hyp==length(opt.pair_indlbls)+1
    overtitle{pp.figdum}='';
    plot_hypchans=1:imp.maxchans;
else
    overtitle{pp.figdum}=sprintf('%s / %s',opt.pair_indlbls{hyp},scl.g_label{group});
    plot_hypinds=find(opt.pair_inds==hyp);
    plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
end
itc_plot_base=permute(meanx(wave_evknormdata(t_start_b:t_end_b,plot_hypchans,:,:,s_inds_g(:,group)),[3 5]),[3 1 2]);
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        itc_plot_data = meanx( bsxfun(@minus, ...
            wave_evknormdata(:,plot_hypchans,:,pp.cond_diff{1},s_inds_g(:,group)), ...
            wave_evknormdata(:,plot_hypchans,:,pp.cond_diff{2},s_inds_g(:,group))), [1 3]);
        v(group,hyp,cond,:) = [minx(itc_plot_data) maxx(itc_plot_data)];
        %[stats, df, pvals, surrog] = statcond( itc_diff_data, 'paired','on','method','perm', ...
            %'naccu', n_perms, 'alpha', p_alpha, 'structoutput', 'on');
        %[p_fdr, p_masked] = fdr( stats.pval, p_alpha, 'nonParametric');
        %[p_masked, p_fdr, adj_ci_cvrg, adj_p] = fdr_bh(pvals, p_alpha,'pdep');
        %itc_plot_data(~stats.mask) = 0;
        %itc_plot_data(~p_masked) = 0;
        %imagesc(fliplr(itc_plot_data)');
        %set(gca,'YDir','normal');
        %pvals_hyp(:,:,hyp)=pvals;
        [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
    else
        itc_diff_data{cond}=bsxfun(@minus,meanx(wave_evknormdata(:,plot_hypchans,:,cond,s_inds_g(:,group)),[1 3 5]), ...
            itc_plot_base);
        %itc_diff_data{cond}=meanx(wave_evknormdata(:,plot_hypchans,:,cond,s_inds_g(:,group)),[1 3 5]);
        itc_plot_data = squeeze(mean( itc_diff_data{cond} , 3));
        [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
        v(group,hyp,cond,:) = caxis;
    end
    %shading flat; 
    %imagesc(ero_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
end
end
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
cmap=makecmap(c, Inf, 'redgreen');
%cmap = parula(256);
%cmap = pmkmp(256, 'cubicl');
cmap_diff=makecmap(c_diff, 0, 'redgreen');
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[1300 120 1500 400]);
for splot=1:length(pp.plotn_cond);
    temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
    if splot==length(pp.plotn_cond)
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
%tightfig;
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle{fig},pp.sp_d);
set_print_size(18,6);
%plottitle(overtitle{fig},2);
end
%make color bars separately
cb_pos=[0.25 0.1 0.05 0.8];
figure;
set(gcf,'position',[2800 120 300 400]);
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
clear_plotassistvars

%% image ITC in time-freq in a REGION, as increase from baseline
% uses condition-mean (common) baseline, GROUPS TOGETHER

sp_rowlabel=scl.g_label{pp.plotn_g};
sp_columnlabel={scl.cond_label{pp.plotn_cond}};
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

spec_chans = [7, 25, 58];
p_alpha=0.05;
n_perms=2000;

%pp.figdum=pp.figdum_init;
pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),6,length(pp.plotn_cond),2);
clear overtitle;
clear itc_stat_data;
for hyp=1:3
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
if hyp==length(opt.pair_indlbls)+1
    overtitle{pp.figdum}='';
    plot_hypchans=1:imp.maxchans;
else
    overtitle{pp.figdum}=sprintf('%s / %s',opt.pair_indlbls{hyp});
    %plot_hypinds=find(opt.pair_inds==hyp);
    %plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
    plot_hypchans = spec_chans(hyp);
    %plot_hypchans = opt.itc_hypchans{hyp};
end

for group=pp.chosen_g(pp.plotn_g)
for cond=pp.plotn_cond(1:end-1)
    itc_stat_data{cond,group-1} = squeeze( bsxfun(@minus,mean(wave_evknormdata(:,plot_hypchans,:,cond,s_inds_g(:,group)),2), ...
            mean(mean(mean( wave_evknormdata(t_start_b:t_end_b,plot_hypchans,:,:,s_inds_g(:,group)), 1),2),4) ) );
end
end
%[stats, df, pvals, surrog] = statcond( itc_stat_data, 'paired', 'off', ...
%        'method', 'perm', 'naccu', 2000, 'alpha', 0.05, 'structoutput', 'on');
%[p_fdr, p_masked] = fdr( stats.pval, p_alpha, 'nonParametric');
%[p_masked, p_fdr, adj_ci_cvrg, adj_p] = fdr_bh(pvals, p_alpha,'pdep');
for group=pp.chosen_g(pp.plotn_g)
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1)*length(pp.plotn_g),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        itc_plot_data = meanx( bsxfun(@minus, ...
            wave_evknormdata(:,plot_hypchans,:,pp.cond_diff{1},s_inds_g(:,group)), ...
            wave_evknormdata(:,plot_hypchans,:,pp.cond_diff{2},s_inds_g(:,group))), [1 3]);
        v(group,hyp,cond,:) = [minx(itc_plot_data) maxx(itc_plot_data)];
        %itc_plot_data(~stats.mask{1}) = 0;
        %itc_plot_data(~p_masked) = 0;
        %imagesc(fliplr(itc_plot_data)');
        %set(gca,'YDir','normal');
        [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
    else
        %itc_stat_data{cond,group-1} = meanx(wave_evknormdata(:,plot_hypchans,:,cond,s_inds_g(:,group)), [1 3]);
        itc_plot_data = squeeze(mean( itc_stat_data{cond,group-1} , 3));
        [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
        v(group,hyp,cond,:) = caxis;
    end
    %shading flat; 
    %imagesc(ero_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label2); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
end
end

c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
cmap=makecmap(c, Inf, 'redgreen');
%cmap = parula(256);
%cmap = pmkmp(256, 'cubicl');
%cmap = pmkmp(256, 'linlhot');
c(1) = 0;
%cmap = flipud(hot(256));
cmap_diff=makecmap(c_diff, 0, 'redgreen');
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[1300 120 1500 600]);
for splot=1:length(pp.plotn_g)*length(pp.plotn_cond);
    temp_s=subplot(pp.sp_d(1)*length(pp.plotn_g),pp.sp_d(2),splot);
    if mod(splot,length(pp.plotn_cond))==0
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
%tightfig;
%adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle{fig},pp.sp_d);
set_print_size(18,12);
%plottitle(overtitle{fig},2);
end
%make color bars separately
figure;
set(gcf,'position',[2800 120 300 400]);
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
%[stats, df, pvals, surrog] = statcond( itc_stat_data, 'paired','off','method','perm', ...
%        'naccu', n_perms, 'alpha', p_alpha, 'structoutput', 'on');

clear_plotassistvars

%% image ITC in time-freq in a REGION, as increase from baseline
% uses condition-mean (common) baseline, GROUPS TOGETHER, GROUP DIFF

sp_rowlabel=scl.g_label(pp.chosen_g(pp.plotn_g));
sp_rowlabel{end+1} = 'Group Contrast';
sp_columnlabel=scl.cond_label(1:imp.maxconds);
%sp_columnlabel=[];
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

spec_chans = [7, 25, 58];
p_alpha=0.05;
n_perms=2000;

%pp.figdum=pp.figdum_init;
pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),6,length(pp.plotn_cond),2);
clear overtitle;
clear itc_stat_data;
for hyp=1
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
if hyp==length(opt.pair_indlbls)+1
    overtitle{pp.figdum}='';
    plot_hypchans=1:imp.maxchans;
else
    overtitle{pp.figdum}=opt.pair_indlbls{hyp};
    %plot_hypinds=find(opt.pair_inds==hyp);
    %plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
    %plot_hypchans = spec_chans(hyp);
    plot_hypchans = opt.itc_hypchans{hyp};
end
gdum=0;
for group=2:4
for cond=1:imp.maxconds
    subplot_dummy=subplot_dummy+1;
    subplot(length(pp.plotn_g)+1,imp.maxconds,subplot_dummy)
    if group==4
        itc_plot_data = meanx(wave_evknormdata(:,plot_hypchans,:,cond,s_inds_g(:,3)),[1 3]) - ...
            meanx(wave_evknormdata(:,plot_hypchans,:,cond,s_inds_g(:,2)),[1 3]);
        % control - alcoholic
        % blue = control has more ITC
        % red = alcoholic has more ITC
        v(group,hyp,cond,:) = [minx(itc_plot_data) maxx(itc_plot_data)];
        [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
    else
        itc_stat_data{cond,group-1} = squeeze( bsxfun(@minus,mean(wave_evknormdata(:,plot_hypchans,:,cond,s_inds_g(:,group)),2), ...
            mean(mean(mean( wave_evknormdata(t_start_b:t_end_b,plot_hypchans,:,:,s_inds_g(:,group)), 1),2),4) ) );
        %itc_stat_data{cond,group-1} = meanx(wave_evknormdata(:,plot_hypchans,:,cond,s_inds_g(:,group)), [1 3]);
        itc_plot_data = squeeze(mean( itc_stat_data{cond,group-1} , 3));
        [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
        set(h,'EdgeColor','None');
        v(group,hyp,cond,:) = caxis;
    end
    %shading flat; 
    %imagesc(ero_plot_data');
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label2); %ylabel('Frequency (Hz)');
    grid on; set(gca,'Layer','Top');
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
end
end
end

c(1)=min(min(min(v(1:end-1,:,:,1)))); c(2)=max(max(max(v(1:end-1,:,:,2))));
c_diff(1)=min(min(v(end,:,:,1))); c_diff(2)=max(max(v(end,:,:,2)));
cmap=makecmap(c, Inf, 'purpleorange');
c(1)=0;
%cmap = parula(256);
%cmap = pmkmp(256, 'cubicl');
cmap_diff=makecmap(c_diff, 0, 'redblue');
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
set(gcf,'position',[1300 120 800 800]);
for splot=1:(length(pp.plotn_g)+1)*imp.maxconds;
    temp_s=subplot(length(pp.plotn_g)+1,imp.maxconds,splot);
    if splot > length(pp.plotn_g)*imp.maxconds
        caxis([c_diff(1) c_diff(2)]);
        colormap(temp_s, cmap_diff);
    else
        caxis([c(1) c(2)]);
        colormap(temp_s, cmap);
    end
end
%tightfig;
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle{fig}, ...
    [length(pp.plotn_g)+1 imp.maxconds]);
set_print_size(12,18);
%plottitle(overtitle{fig},2);
end
%make color bars separately
cb_pos=[0.25 0.1 0.05 0.8];
figure;
set(gcf,'position',[2800 120 300 400]);
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
%[stats, df, pvals, surrog] = statcond( itc_stat_data, 'paired','off','method','perm', ...
%        'naccu', n_perms, 'alpha', p_alpha, 'structoutput', 'on');

clear_plotassistvars

%% error bar plot of a TFROI for ITC in a region
% original

t_ms = [200 400];
f_hz = [2 3];

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

[~,t_start] = min( abs(scl.t_ms - t_ms(1) ));
[~,t_end] = min( abs(scl.t_ms - t_ms(2) ));

[~,f_start] = min( abs(scl.freqs - f_hz(1) ));
[~,f_end] = min( abs(scl.freqs - f_hz(2) ));

bar_plot_y = zeros(1,2);
bar_plot_e = zeros(1,2);

for hyp=1:3
%figure;
plot_hypinds=find(opt.pair_inds==hyp);
plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));

for group=pp.plotn_g

itc_plot_base = meanx( wave_evknormdata(t_start_b:t_end_b, plot_hypchans, ...
    f_end:f_start, :, s_inds_g(:,group)), 5 );
clear itc_diff_data
for cond=1:2
    itc_diff_data{cond}=bsxfun(@minus, ...
        meanx( wave_evknormdata(t_start:t_end, plot_hypchans, f_end:f_start, ...
            cond, s_inds_g(:,group)), 5 ), ...
        itc_plot_base);
    bar_plot_y(cond,hyp) = mean( itc_diff_data{cond} );
    bar_plot_e(cond,hyp) = std( itc_diff_data{cond} ) / sqrt(sum( s_inds_g(:,group) ));
end
%errorbarbar(bar_plot_y, bar_plot_e);
%bar(bar_plot_y(hyp,:)); hold on;
%errorb(bar_plot_y(hyp,:), bar_plot_e(hyp,:));
%hold off;
end
end
%set_print_size
%tightfig;
clear_plotassistvars

%% error bar plot of a TFROI for ITC in a region
% groups

t_ms = [200 400];
f_hz = [4 6];

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

[~,t_start] = min( abs(scl.t_ms - t_ms(1) ));
[~,t_end] = min( abs(scl.t_ms - t_ms(2) ));

[~,f_start] = min( abs(scl.freqs - f_hz(1) ));
[~,f_end] = min( abs(scl.freqs - f_hz(2) ));

bar_plot_y = zeros(1,2);
bar_plot_e = zeros(1,2);

plotdata_cell = cell(2,1);
for freq_range=1 %:3
[~,f_start] = min( abs(scl.freqs - f_hz(freq_range) ) );
[~,f_end] = min( abs(scl.freqs - f_hz(freq_range) ) );
for hyp=1:3
plot_hypinds=find(opt.pair_inds==hyp);
plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
for group=pp.chosen_g(pp.plotn_g)
    
    % ITC regular
    %plotdata_cell{group-1, freq_range} = [ plotdata_cell{group-1, freq_range}, ...
    %    meanx(wave_evknormdata(t_start:t_end,plot_hypchans,f_end:f_start,:,s_inds_g(:,group)), [4 5])' ];
    
    % ITC baselined
    plotdata_cell{group-1, freq_range} = [ plotdata_cell{group-1, freq_range}, meanx(bsxfun(@minus, ...
        wave_evknormdata(t_start:t_end,plot_hypchans,f_end:f_start,:,s_inds_g(:,group)), ...
        mean(mean(wave_evknormdata(t_start_b:t_end_b,plot_hypchans,f_end:f_start,:,s_inds_g(:,group)),1),4)),[4 5])' ];
    
    % ICPS regular
    %plotdata_cell{group-1, freq_range} = [ plotdata_cell{group-1, freq_range}, ...
    %    meanx(cohdata(t_start:t_end,f_end:f_start,:,plot_hypinds,s_inds_g(:,group)), [3 5])' ];
    
    % ICPS baselined
    %plotdata_cell{group-1, freq_range} = [ plotdata_cell{group-1, freq_range}, meanx(bsxfun(@minus, ...
    %    wave_evknormdata(t_start:t_end,f_end:f_start,:,plot_hypinds,s_inds_g(:,group)), ...
    %    mean(mean(wave_evknormdata(t_start_b:t_end_b,f_end:f_start,:,plot_hypinds,s_inds_g(:,group)),1),3) ) ,[3 5])' ];
    
    %bar_plot_y(cond,group,hyp) = mean( itc_stat_data{cond,group-1} );
    %bar_plot_e(cond,group,hyp) = std( itc_stat_data{cond,group-1} ) / sqrt(sum( s_inds_g(:,group) ));

end
end
end
clear_plotassistvars

%% distribution of p-values

figure;
subplot_dummy=0;
for hyp=[1 2 3 7]
    subplot_dummy=subplot_dummy+1;
    subplot(2,2,subplot_dummy);
    %hist( reshape( pvals_hyp(:,:,hyp), [1 prod(size(pvals))] ) )
    [f,xi] = ksdensity( reshape( pvals_hyp(:,:,hyp), [1 prod(size(pvals))] ) );
    plot(xi,f)
    axis([0 1 0 3])
    %ylim([0 700])
end


%% image ITC in time-freq at a chosen channel, with ERP superimposed

sp_rowlabel={''};
sp_columnlabel=scl.cond_label;
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';
subplot_dims=pp.sp_d;

pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
for group=pp.chosen_g
for chan=pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
overtitle{pp.figdum}=sprintf('%s / %s',scl.chan_label{chan},scl.g_label{group});
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
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
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
    grid on;
    %title(['ITC at ',scl.chan_label{chan},' : ',scl.g_label{group},' : ',scl.cond_label{cond}])
    hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--');
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
end
end
c(1)=min(min(min(v(:,:,1:end-1,1)))); c(2)=max(max(max(v(:,:,1:end-1,2))));
c_diff(1)=min(min(v(:,:,end,1))); c_diff(2)=max(max(v(:,:,end,2)));
for fig=pp.figdum_init+1:pp.figdum
figure(fig)
for splot=1:length(pp.plotn_cond);
    if splot==length(pp.plotn_cond)
        subplot(pp.sp_d(1),pp.sp_d(2),splot); caxis([c_diff(1) c_diff(2)]);
    else
        subplot(pp.sp_d(1),pp.sp_d(2),splot); caxis([c(1) c(2)]);
    end
    colorbar;
end
plottitle(overtitle{fig});
tightfig;
set_print_size(20,8);
end
%distFig('s','ext','transpose',true);
clear_plotassistvars

%% plot ITC as a topographic plot over the head (after subtracting common baseline)

%warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)-1];

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

itc_topo_scale=[0 .37];
itc_diff_limits=[-.05 .1];
cmap=makecmap(itc_topo_scale,Inf,'redgreen',256);
cmap_diff=makecmap(itc_diff_limits,0,'redgreen',256);

v=zeros(length(pp.plotn_cond),length(pp.chosen_g),length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms),2);
clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for freq_range=2
figure; subplot_dummy=0;
overtitle=sprintf('Topography of ITC in %1.1f - %1.1f Hz / %s',pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.g_label{group});
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
itc_topo_base=meanx(wave_evknormdata(t_start_b:t_end_b,:,f_end:f_start,:,s_inds_g(:,group)), 2);
%itc_topo_base=zeros(size(meanx(wave_evknormdata(t_start_b:t_end_b,:,f_end:f_start,:,s_inds_g(:,group)), 2)));
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)-1
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        subplot_dummy=subplot_dummy+1;
        sp_temp=subplot(length(pp.plotn_cond),length(pp.t_start_ms)-1,subplot_dummy);
        if cond==imp.maxconds+1
            topo_data=meanx(wave_evknormdata(t_start:t_end,:,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),2) - ...
                meanx(wave_evknormdata(t_start:t_end,:,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),2);
            h=topoplot(topo_data,chan_locs,'maplimits',[itc_diff_limits(1) itc_diff_limits(2)], ...
                'electrodes',pp.topo_elecs,'colormap',cmap_diff,'style','both','numcontour',7);
            colormap(sp_temp,cmap_diff);
        else
            topo_data=meanx(wave_evknormdata(t_start:t_end,:,f_end:f_start,cond,s_inds_g(:,group)),2) - ...
                itc_topo_base;
            h=topoplot(topo_data,chan_locs,'maplimits',[itc_topo_scale(1) itc_topo_scale(2)], ...
                'electrodes',pp.topo_elecs,'colormap',cmap,'style','both','numcontour',7);
            colormap(sp_temp,cmap);
        end
        %set(h,'EdgeColor','None');
        v(cond,group,freq_range,win,1) = min(topo_data); v(cond,group,freq_range,win,2) = max(topo_data);        
    end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims,false);
%tightfig;
set_print_size(12,12);
end
end
c(1)=min(min(min(min(v(1:end-1,:,:,:,1))))); c(2)=max(max(max(max(v(1:end-1,:,:,:,2)))));
c_diff(1)=min(min(min(v(end,:,:,:,1)))); c_diff(2)=max(max(max(v(end,:,:,:,2))));
figure;
colorscale_plot(itc_topo_scale, cmap, 0.25);
colorscale_plot(itc_diff_limits, cmap_diff, 0.75);
clear_plotassistvars;

%% image ITC as a scalp plot in a series of time-frequency windows

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.plotn_cond),length(pp.t_start_ms)];

itc_topo_scale=[0.2 0.5];
itc_diff_limits=[-.03 .08];
cmap=makecmap(itc_topo_scale);
cmap_diff=makecmap(itc_diff_limits);

%
v=zeros(length(pp.plotn_cond),length(pp.chosen_g),length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms),2);
for group=pp.chosen_g(pp.plotn_g)
for freq_range=2 %pp.plotn_f
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
        if win~=pp.maxwin
        if cond==imp.maxconds+1
            topo_data=squeeze(mean(mean(mean(mean(itcdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4),5)) - ...
                squeeze(mean(mean(mean(mean(itcdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4),5));
            topoplot(topo_data,chan_locs,'maplimits',[itc_diff_limits(1) itc_diff_limits(2)],...
                'electrodes',pp.topo_elecs,'colormap',cmap_diff,'style','fill','numcontour',7);
            freezeColors;
        else
            topo_data=squeeze(mean(mean(mean(itcdata(t_start:t_end,pp.chosen_topochan,f_end:f_start,cond,s_inds_g(:,group)),1),3),5));
            topoplot(topo_data,chan_locs,'maplimits',[itc_topo_scale(1) itc_topo_scale(2)],...
                'electrodes',pp.topo_elecs,'colormap',cmap,'style','fill','numcontour',7);
            freezeColors;
        end
        shading flat
        else
            set(gca,'Visible','Off');
        end
        v(cond,group,freq_range,win,1) = min(topo_data); v(cond,group,freq_range,win,2) = max(topo_data);
        %title(sprintf('%s, %d - %d ms, %1.1f - %1.1f Hz, %s',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),pp.f_start_hz(freq_range),pp.f_end_hz(freq_range),scl.g_label{group}))
    end
    %color bar
    cb_pos = get(gca,'Position') + [0.08 0 0 0];
    cb_pos(3) = 0.05; %set width
    if cond==imp.maxconds+1
        cb_lim = itc_diff_limits;
    else
        cb_lim = itc_topo_scale;
    end
    colorscale([1 256], cb_lim, range(cb_lim)/5, 'vert', ...
        'Position',cb_pos);
    freezeColors;
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
%fig(pp.figdum, 'units', 'centimeters', 'width', 17, 'height', length(pp.plotn_cond)*3, ...
%    'font', 'Monospaced', 'fontsize', 11);
set_print_size(17,17/scl.phi);
end
end
c(1)=min(min(min(min(v(1:end-1,:,:,:,1))))); c(2)=max(max(max(max(v(1:end-1,:,:,:,2)))));
c_diff(1)=min(min(min(v(end,:,:,:,1)))); c_diff(2)=max(max(max(v(end,:,:,:,2))));
clear_plotassistvars

%% ITC - create condition bar plots with error bars

sp_rowlabel=make_freqlabels(pp.f_start_hz(pp.plotn_f),pp.f_end_hz(pp.plotn_f));
sp_columnlabel=make_timelabels(pp.t_start_ms,pp.t_end_ms);
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[length(pp.f_start_hz(pp.plotn_f)),length(pp.t_start_ms)];

width=[];
bw_xlabel=[];
bw_ylabel=[];
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
        axis([0 4 0.15 0.6])
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

%% plot per-subject lines of the shape of results for a certain measure across conditions (ERO)

chosen_chan=7;
%win_t=[200 400];
win_t=[300 500];
win_f=[4 5.3];

colorkey={'r','g','b','m','k','c'};
sitekey={'UConn','Indiana','Iowa','SUNY','WashU','UCSD'};
counter=zeros(1,6);

[~,t_s]=min(abs(scl.t_ms-win_t(1))); [~,t_e]=min(abs(scl.t_ms-win_t(2)));
[~,f_s]=min(abs(scl.freqs-win_f(1))); [~,f_e]=min(abs(scl.freqs-win_f(2)));

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

figure; gdum=0;
overtitle=sprintf('ERSP at %s in %1.1f - %1.1f Hz from %d - %d ms',...
    scl.chan_label{chosen_chan},win_f(1),win_f(2),win_t(1),win_t(2));
for group=pp.chosen_g(pp.plotn_g)
    conddiff_dir=zeros(imp.maxconds,1);
    gdum=gdum+1;
    %subplot(1,2,gdum)
    linedata = meanx( bsxfun (@rdivide, wave_totpowdata(t_s:t_e,chosen_chan,f_e:f_s,:,s_inds_g(:,group)), ...
        mean(wave_totpowdata(t_start_b:t_end_b,chosen_chan,f_e:f_s,:,s_inds_g(:,group)),1) ) , [4 5] );
    %linedata=meanx(wave_totpowdata(t_s:t_e,chosen_chan,f_e:f_s,:,s_inds_g(:,group)),[4 5]) - ...
    %    meanx(wave_totpowdata(1:scl.t_zero,chosen_chan,f_e:f_s,:,s_inds_g(:,group)),[4 5]); %ERSP
    ranova_data{gdum}=linedata';
    for s=1:size(linedata,2)
    %site=str2double(mat_list{s}(65));
    %counter(site)=counter(site)+1;
    if linedata(pp.cond_diff{1},s) > linedata(pp.cond_diff{2},s)
        mark=['b-o'];
        conddiff_dir(1)=conddiff_dir(1)+1;
    else
        mark=['m--*'];
        conddiff_dir(2)=conddiff_dir(2)+1;
    end
    %subplot(1,6,site);
    %title(sitekey{site});
    plot(linedata(:,s),mark); hold on; axis([0 imp.maxconds+1 0.5 2.5]);
    end
    conddiff_text=sprintf(['%s > %s : %d/%d (%2.0f%%), M=%2.1f',char(177),'%1.1f'],...
        scl.cond_label{pp.cond_diff{1}},scl.cond_label{pp.cond_diff{2}},...
        conddiff_dir(1),sum(s_inds_g(:,group)),...
        conddiff_dir(1)/sum(s_inds_g(:,group))*100,...
        abs(diff(mean(linedata(:,linedata(pp.cond_diff{1},:)>...
        linedata(pp.cond_diff{2},:)),2))),...
        abs(diff(std(linedata(:,linedata(pp.cond_diff{1},:)>...
        linedata(pp.cond_diff{2},:)),0,2))));
    text(0.2,2.2,conddiff_text)
    set(gca,'XTick',[1:imp.maxconds],'XTickLabel',scl.cond_label);
    xlabel('Condition');
    ylabel('ERSP');
    %title(scl.g_label{group});    
end
[p,ranova_table]=anova_rm(ranova_data,'off');
plottitle(overtitle);
tightfig;
clear_plotassistvars;

%% plot per-subject lines of the shape of results for a certain measure across conditions (ITC)

chosen_chan=7;
%win_t=[200 400];
win_t=[200 400];
win_f=[4 6];

colorkey={'r','g','b','m','k','c'};
sitekey={'UConn','Indiana','Iowa','SUNY','WashU','UCSD'};
counter=zeros(1,6);

[~,t_s]=min(abs(scl.t_ms-win_t(1))); [~,t_e]=min(abs(scl.t_ms-win_t(2)));
[~,f_s]=min(abs(scl.freqs-win_f(1))); [~,f_e]=min(abs(scl.freqs-win_f(2)));

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

gdum=0;
overtitle=sprintf('ITC at %s in %1.1f - %1.1f Hz from %d - %d ms',...
    scl.chan_label{chosen_chan},win_f(1),win_f(2),win_t(1),win_t(2));
for group=3 %pp.chosen_g(pp.plotn_g)
    figure;
    conddiff_dir=zeros(imp.maxconds,1);
    gdum=gdum+1;
    %subplot(1,2,gdum)
    linedata = meanx( bsxfun (@minus, wave_evknormdata(t_s:t_e,chosen_chan,f_e:f_s,:,s_inds_g(:,group)), ...
        mean(wave_evknormdata(t_start_b:t_end_b,chosen_chan,f_e:f_s,:,s_inds_g(:,group)),1) ) , [4 5] );
    %linedata=meanx(wave_totpowdata(t_s:t_e,chosen_chan,f_e:f_s,:,s_inds_g(:,group)),[4 5]) - ...
    %    meanx(wave_totpowdata(1:scl.t_zero,chosen_chan,f_e:f_s,:,s_inds_g(:,group)),[4 5]); %ERSP
    ranova_data{gdum}=linedata';
    for s=1:size(linedata,2)
    %site=str2double(mat_list{s}(65));
    %counter(site)=counter(site)+1;
    if linedata(pp.cond_diff{1},s) > linedata(pp.cond_diff{2},s)
        mark=['b-o'];
        conddiff_dir(1)=conddiff_dir(1)+1;
    else
        mark=['m--*'];
        conddiff_dir(2)=conddiff_dir(2)+1;
    end
    %subplot(1,6,site);
    %title(sitekey{site});
    plot(linedata(:,s),mark); hold on; axis([0 imp.maxconds+1 -.2 .6]);
    end
    conddiff_text=sprintf(['%s > %s : %d/%d (%2.0f%%), M=%2.1f',char(177),'%1.2f'],...
        scl.cond_label{pp.cond_diff{1}},scl.cond_label{pp.cond_diff{2}},...
        conddiff_dir(1),sum(s_inds_g(:,group)),...
        conddiff_dir(1)/sum(s_inds_g(:,group))*100,...
        abs(diff(mean(linedata(:,linedata(pp.cond_diff{1},:)>...
        linedata(pp.cond_diff{2},:)),2))),...
        abs(diff(std(linedata(:,linedata(pp.cond_diff{1},:)>...
        linedata(pp.cond_diff{2},:)),0,2))));
    text(0.2,0.5,conddiff_text)
    set(gca,'XTick',[1:imp.maxconds],'XTickLabel',scl.cond_label);
    xlabel('Condition');
    ylabel('ITC');
    %title(scl.g_label{group});    
end
[p,ranova_table]=anova_rm(ranova_data,'off');
plottitle(overtitle);
tightfig;

diff_bysub = diff(linedata);

clear_plotassistvars;

%% plot per-subject lines of the shape of results for a certain measure across conditions (ITC)
% region, 1 > 2

hyp=1;
plot_hypinds=opt.pair_inds==hyp;
plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
%win_t=[200 400];
win_t=[200 400];
win_f=[4 7];

[~,t_s]=min(abs(scl.t_ms-win_t(1))); [~,t_e]=min(abs(scl.t_ms-win_t(2)));
[~,f_s]=min(abs(scl.freqs-win_f(1))); [~,f_e]=min(abs(scl.freqs-win_f(2)));

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

gdum=0;
overtitle=sprintf('ITC at %s in %1.1f - %1.1f Hz from %d - %d ms',...
    opt.pair_indlbls{hyp},win_f(1),win_f(2),win_t(1),win_t(2));
for group=pp.chosen_g(pp.plotn_g)
    figure;
    conddiff_dir=zeros(imp.maxconds,1);
    gdum=gdum+1;
    %subplot(1,2,gdum)
    linedata = meanx( bsxfun (@minus, wave_evknormdata(t_s:t_e,plot_hypchans,f_e:f_s,:,s_inds_g(:,group)), ...
        mean(mean(wave_evknormdata(t_start_b:t_end_b,plot_hypchans,f_e:f_s,:,s_inds_g(:,group)),1),4) ) , [4 5] );
    ranova_data{gdum}=linedata';
    for s=1:size(linedata,2)
    if linedata(pp.cond_diff{1},s) > linedata(pp.cond_diff{2},s)
        mark=['r-']; sp=1;
        conddiff_dir(1)=conddiff_dir(1)+1;
    else
        mark=['g-']; sp=2;
        conddiff_dir(2)=conddiff_dir(2)+1;
    end
    subplot(1,2,sp);
    plot(linedata(:,s),mark); hold on; axis([0.75 imp.maxconds+.25 -.2 .6]);
    end
    conddiff_text=sprintf(['%s > %s : %d/%d (%2.0f%%), M=%2.1f',char(177),'%1.2f'],...
        scl.cond_label{pp.cond_diff{1}},scl.cond_label{pp.cond_diff{2}},...
        conddiff_dir(1),sum(s_inds_g(:,group)),...
        conddiff_dir(1)/sum(s_inds_g(:,group))*100,...
        abs(diff(mean(linedata(:,linedata(pp.cond_diff{1},:)>...
        linedata(pp.cond_diff{2},:)),2))),...
        abs(diff(std(linedata(:,linedata(pp.cond_diff{1},:)>...
        linedata(pp.cond_diff{2},:)),0,2))));
    text(1,0.5,conddiff_text)
    set(gca,'XTick',[1:imp.maxconds],'XTickLabel',scl.cond_label);
    xlabel('Condition');
    ylabel('ITC');
    title(scl.g_label{group});    
end
[p,ranova_table]=anova_rm(ranova_data,'off');
plottitle(overtitle);
set_print_size(4,6);

diff_bysub = diff(linedata);

%clear_plotassistvars;

%% plot per-subject lines of the shape of results for a certain measure across conditions (ITC)
% region, 2 > 1

hyp=2;
plot_hypinds=find(opt.pair_inds==hyp);
plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
%win_t=[200 400];
win_t=[200 400];
win_f=[2 3];

[~,t_s]=min(abs(scl.t_ms-win_t(1))); [~,t_e]=min(abs(scl.t_ms-win_t(2)));
[~,f_s]=min(abs(scl.freqs-win_f(1))); [~,f_e]=min(abs(scl.freqs-win_f(2)));

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

figure; gdum=0;
overtitle=sprintf('ITC at %s in %1.1f - %1.1f Hz from %d - %d ms',...
    opt.pair_indlbls{hyp},win_f(1),win_f(2),win_t(1),win_t(2));
for group=pp.chosen_g(pp.plotn_g)
    conddiff_dir=zeros(imp.maxconds,1);
    gdum=gdum+1;
    %subplot(1,2,gdum)
    linedata = meanx( bsxfun (@minus, wave_evknormdata(t_s:t_e,plot_hypchans,f_e:f_s,:,s_inds_g(:,group)), ...
        mean(mean(wave_evknormdata(t_start_b:t_end_b,plot_hypchans,f_e:f_s,:,s_inds_g(:,group)),1),4) ) , [4 5] );
    ranova_data{gdum}=linedata';
    for s=1:size(linedata,2)
    if linedata(pp.cond_diff{1},s) < linedata(pp.cond_diff{2},s)
        mark=['b-o'];
        conddiff_dir(1)=conddiff_dir(1)+1;
    else
        mark=['m--*'];
        conddiff_dir(2)=conddiff_dir(2)+1;
    end
    plot(linedata(:,s),mark); hold on; axis([0.75 imp.maxconds+.25 -.2 .6]);
    end
    conddiff_text=sprintf(['%s > %s : %d/%d (%2.0f%%), M=%2.1f',char(177),'%1.2f'],...
        scl.cond_label{pp.cond_diff{2}},scl.cond_label{pp.cond_diff{1}},...
        conddiff_dir(1),sum(s_inds_g(:,group)),...
        conddiff_dir(1)/sum(s_inds_g(:,group))*100,...
        abs(diff(mean(linedata(:,linedata(pp.cond_diff{1},:)<...
        linedata(pp.cond_diff{2},:)),2))),...
        abs(diff(std(linedata(:,linedata(pp.cond_diff{1},:)<...
        linedata(pp.cond_diff{2},:)),0,2))));
    text(1,0.5,conddiff_text)
    set(gca,'XTick',[1:imp.maxconds],'XTickLabel',scl.cond_label);
    xlabel('Condition');
    ylabel('ITC');
    %title(scl.g_label{group});
end
[p,ranova_table]=anova_rm(ranova_data,'off');
plottitle(overtitle);

diff_bysub = diff(linedata);

clear_plotassistvars;

%% boxplot conditions

hyp=2;
plot_hypinds=opt.pair_inds==hyp;
plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
%win_t=[200 400];
win_t=[200 400];
win_f=[2 3];

[~,t_s]=min(abs(scl.t_ms-win_t(1))); [~,t_e]=min(abs(scl.t_ms-win_t(2)));
[~,f_s]=min(abs(scl.freqs-win_f(1))); [~,f_e]=min(abs(scl.freqs-win_f(2)));

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

figure; gdum=0;
overtitle=sprintf('ITC at %s in %1.1f - %1.1f Hz from %d - %d ms',...
    opt.pair_indlbls{hyp},win_f(1),win_f(2),win_t(1),win_t(2));
for group=pp.chosen_g(pp.plotn_g)
    conddiff_dir=zeros(imp.maxconds,1);
    gdum=gdum+1;
    %subplot(1,2,gdum)
    linedata = meanx( bsxfun (@minus, wave_evknormdata(t_s:t_e,plot_hypchans,f_e:f_s,:,s_inds_g(:,group)), ...
        mean(mean(wave_evknormdata(t_start_b:t_end_b,plot_hypchans,f_e:f_s,:,s_inds_g(:,group)),1),4) ) , [4 5] );
    ranova_data{gdum}=linedata';
    boxplot(linedata');
    ylim([-.05 .4]);
    text(0.75,0.35,conddiff_text)
    set(gca,'XTick',[1:imp.maxconds],'XTickLabel',scl.cond_label);
    xlabel('Condition');
    ylabel('ITC');
    %title(scl.g_label{group});    
end
[p,ranova_table]=anova_rm(ranova_data,'off');
plottitle(overtitle);
set_print_size(4,6);

diff_bysub = diff(linedata);

clear_plotassistvars;

%% scatter baseline with event-related ITC

win_t=[400 600];

[~,t_s]=min(abs(scl.t_ms-win_t(1)));
[~,t_e]=min(abs(scl.t_ms-win_t(2)));

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

for group=pp.chosen_g(pp.plotn_g)
    for chan=7
        figure;
        subplot_dummy=0;
        for freq_range=1:3
            subplot(1,3,freq_range);
            [~,f_s]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
            [~,f_e]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
            
            itc_plot_base = meanx( wave_evknormdata(t_start_b:t_end_b,chan,f_e:f_s,:,s_inds_g(:,group)), 5);
            itc_plot_data = meanx( wave_evknormdata(t_s:t_e,chan,f_e:f_s,:,s_inds_g(:,group)), 5);
            scatter(itc_plot_base, itc_plot_data); hold on;
            plot(linspace(0,1,100),linspace(0,1,100),'k--'); hold off;
        end
    end
end

clear_plotassistvars;

%% plot per-subject lines of the shape of results for a certain measure across conditions (ITC)
% region, divisive to the first condition ???

hyp=1;
plot_hypinds=find(opt.pair_inds==hyp);
plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
%win_t=[200 400];
win_t=[200 400];
win_f=[4 6];

[~,t_s]=min(abs(scl.t_ms-win_t(1))); [~,t_e]=min(abs(scl.t_ms-win_t(2)));
[~,f_s]=min(abs(scl.freqs-win_f(1))); [~,f_e]=min(abs(scl.freqs-win_f(2)));

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

figure; gdum=0;
overtitle=sprintf('ITC at %s in %1.1f - %1.1f Hz from %d - %d ms',...
    opt.pair_indlbls{hyp},win_f(1),win_f(2),win_t(1),win_t(2));
for group=pp.chosen_g(pp.plotn_g)
    conddiff_dir=zeros(imp.maxconds,1);
    gdum=gdum+1;
    %subplot(1,2,gdum)
    linedata = meanx( bsxfun (@minus, wave_evknormdata(t_s:t_e,plot_hypchans,f_e:f_s,:,s_inds_g(:,group)), ...
        mean(wave_evknormdata(t_start_b:t_end_b,plot_hypchans,f_e:f_s,:,s_inds_g(:,group)),1) ) , [4 5] );
    %linedata=meanx(wave_totpowdata(t_s:t_e,chosen_chan,f_e:f_s,:,s_inds_g(:,group)),[4 5]) - ...
    %    meanx(wave_totpowdata(1:scl.t_zero,chosen_chan,f_e:f_s,:,s_inds_g(:,group)),[4 5]); %ERSP
    ranova_data{gdum}=linedata';
    for s=1:size(linedata,2)
    %site=str2double(mat_list{s}(65));
    %counter(site)=counter(site)+1;
    if linedata(pp.cond_diff{1},s) > linedata(pp.cond_diff{2},s)
        mark=['b-o'];
        conddiff_dir(1)=conddiff_dir(1)+1;
    else
        mark=['m--*'];
        conddiff_dir(2)=conddiff_dir(2)+1;
    end
    %subplot(1,6,site);
    %title(sitekey{site});
    plot(linedata(:,s),mark); hold on; axis([0 imp.maxconds+1 -.2 .6]);
    end
    conddiff_text=sprintf(['%s > %s : %d/%d (%2.0f%%), M=%2.1f',char(177),'%1.2f'],...
        scl.cond_label{pp.cond_diff{1}},scl.cond_label{pp.cond_diff{2}},...
        conddiff_dir(1),sum(s_inds_g(:,group)),...
        conddiff_dir(1)/sum(s_inds_g(:,group))*100,...
        abs(diff(mean(linedata(:,linedata(pp.cond_diff{1},:)>...
        linedata(pp.cond_diff{2},:)),2))),...
        abs(diff(std(linedata(:,linedata(pp.cond_diff{1},:)>...
        linedata(pp.cond_diff{2},:)),0,2))));
    text(0.2,0.5,conddiff_text)
    set(gca,'XTick',[1:imp.maxconds],'XTickLabel',scl.cond_label);
    xlabel('Condition');
    ylabel('ITC');
    %title(scl.g_label{group});    
end
[p,ranova_table]=anova_rm(ranova_data,'off');
plottitle(overtitle);
tightfig;

diff_bysub = diff(linedata);

clear_plotassistvars;

%% scatter per-subject ITC with behavioral criteria

chosen_chan=7;
win_t=[200 400];
win_f=[4 5.3];

[~,t_s]=min(abs(scl.t_ms-win_t(1))); [~,t_e]=min(abs(scl.t_ms-win_t(2)));
[~,f_s]=min(abs(scl.freqs-win_f(1))); [~,f_e]=min(abs(scl.freqs-win_f(2)));

fitmat=zeros(2,4);

for po=1:4

figure;
overtitle=sprintf('ITC at %s in %1.1f - %1.1f Hz from %d - %d ms',...
    scl.chan_label{chosen_chan},win_f(1),win_f(2),win_t(1),win_t(2));

%xdata=behdata.crit(s_inds_g(:,scl.g_all));
%xdata=behdata.avgbet(s_inds_g(:,scl.g_all));
xdata=behdata.avgbet_po(po,s_inds_g(:,scl.g_all));
%xdata=behdata.stdbet_po(po,s_inds_g(:,scl.g_all));
%xdata=behdata.avgbet_po2(po,s_inds_g(:,scl.g_all));
ydata=meanx(itcdata(t_s:t_e,chosen_chan,f_e:f_s,pp.cond_diff{1},s_inds_g(:,scl.g_all)),5) - ...
    meanx(itcdata(t_s:t_e,chosen_chan,f_e:f_s,pp.cond_diff{2},s_inds_g(:,scl.g_all)),5);
%ydata=meanx(itcdata(t_s:t_e,chosen_chan,f_e:f_s,pp.cond_diff{1},s_inds_g(:,scl.g_all)),5);

fitmat(:,po)=robustfit(xdata,ydata);

scatter(xdata,ydata); hold on;
plot(linspace(10,50,100),linspace(10,50,100)*fitmat(2,po)+fitmat(1,po))


%set(gca,'XTick',[1:imp.maxconds],'XTickLabel',scl.cond_label);
xlabel('Criterion');
ylabel('ITC');
%title(scl.g_label{group});    

plottitle(overtitle);
end
clear_plotassistvars

%% scatter per-subject ERO with behavioral criteria

chosen_chan=7;
win_t=[250 500];
win_f=[4.6 6.4];

[~,t_s]=min(abs(scl.t_ms-win_t(1))); [~,t_e]=min(abs(scl.t_ms-win_t(2)));
[~,f_s]=min(abs(scl.freqs-win_f(1))); [~,f_e]=min(abs(scl.freqs-win_f(2)));

for po=1:4

figure;
overtitle=sprintf('ITC at %s in %1.1f - %1.1f Hz from %d - %d ms',...
    scl.chan_label{chosen_chan},win_f(1),win_f(2),win_t(1),win_t(2));

for group=pp.chosen_g

%xdata=behdata.crit(s_inds_g(:,scl.g_all));
%xdata=behdata.avgbet(s_inds_g(:,scl.g_all));
xdata=behdata.crit_po(po,s_inds_g(:,group));
%xdata=behdata.avgbet_po(po,s_inds_g(:,scl.g_all));
%data=behdata.avgbet_po2(po,s_inds_g(:,scl.g_all));
ydata=meanx(wave_totpowdata(t_s:t_e,chosen_chan,f_e:f_s,pp.cond_diff{1},s_inds_g(:,group)),5) - ...
    meanx(wave_totpowdata(t_s:t_e,chosen_chan,f_e:f_s,pp.cond_diff{2},s_inds_g(:,group)),5);
%ydata=meanx(wave_totpowdata(t_s:t_e,chosen_chan,f_e:f_s,pp.cond_diff{1},s_inds_g(:,group)),5);

scatter(xdata,ydata,scl.g_color{group}); hold on;

end

%set(gca,'XTick',[1:imp.maxconds],'XTickLabel',scl.cond_label);
xlabel('Criterion');
ylabel('ERO');
%title(scl.g_label{group});    

plottitle(overtitle);
end
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
    plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(1),itc_axes(2),100)*...
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
                plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(1),itc_axes(2),100)*p_tfwin(2,chan,freq_range,win)+p_tfwin(1,chan,freq_range,win), scl.g_color{group}); hold on;
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

ero_axes=[min(s_demogs.age_eeg) max(s_demogs.age_eeg) 10 100];
ero_axes_diff=[min(s_demogs.age_eeg) max(s_demogs.age_eeg) -30 20];

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_start_hz(pp.plotn_f)), ...
    length(pp.plotn_cond),pp.maxwin,length(pp.plotn_g));
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
        ero_scatterdata_y=squeeze(mean(mean(mean(wave_totpowdata(...
            t_start:t_end,chan,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),1),3),4)) -...
            squeeze(mean(mean(mean(wave_totpowdata(...
            t_start:t_end,chan,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),1),3),4));
        axis(ero_axes_diff);
    else
        ero_scatterdata_y=squeeze(mean(mean(mean(wave_totpowdata(...
            t_start:t_end,chan,f_end:f_start,cond,s_inds_g(:,group)),1),3),4));
        axis(ero_axes);
    end

    p_tfwin(:,chan,freq_range,cond,win,group)=robustfit(ero_scatterdata_x,ero_scatterdata_y);
    %
    scatter_h(group)=scatter(ero_scatterdata_x,ero_scatterdata_y, scl.g_color{group}); hold on;
    %
    %plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100),'k--'); hold on;
    plot(linspace(ero_axes(1),ero_axes(2),100),linspace(ero_axes(1),ero_axes(2),100)*...
        p_tfwin(2,chan,freq_range,cond,win,group)+p_tfwin(1,chan,freq_range,cond,win,group), scl.g_color{group}); hold on;
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

itc_axes=[min(s_demogs.age_eeg) max(s_demogs.age_eeg) 0 1];
itc_axes_diff=[min(s_demogs.age_eeg) max(s_demogs.age_eeg) -0.5 0.5];

p_tfwin=zeros(2,length(pp.chosen_chan(pp.plotn_chan)),length(pp.f_start_hz(pp.plotn_f)), ...
    length(pp.plotn_cond),pp.maxwin,length(pp.plotn_g));
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

    %p_tfwin(:,chan,freq_range,cond,win,group)=robustfit(itc_scatterdata_x,itc_scatterdata_y);
    p_tfwin(:,chan,freq_range,cond,win,group)=polyfit(itc_scatterdata_x,itc_scatterdata_y,1);
    %
    scatter_h(group)=scatter(itc_scatterdata_x,itc_scatterdata_y, scl.g_color{group}); hold on;
    %
    %plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(3),itc_axes(4),100),'k--'); hold on;
    plot(linspace(itc_axes(1),itc_axes(2),100),linspace(itc_axes(1),itc_axes(2),100)*...
        p_tfwin(1,chan,freq_range,cond,win,group)+p_tfwin(2,chan,freq_range,cond,win,group), scl.g_color{group}); hold on;
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
        contourf(itc_freqcorr,pp.n_contour); %caxis(rho_limits);
        colormap(pp.cmap)
        v(cond,win,:)=caxis;
        title(sprintf('%s, %d - %d ms, %s',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win),scl.chan_label{chan}));
        %title(sprintf('%s, %d - %d ms',scl.cond_label{cond},pp.t_start_ms(win),pp.t_end_ms(win)));
        xlabel('Freq. (Hz)'); ylabel('Freq. (Hz)');
        set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label);
        set(gca,'XTick',scl.f_ytick,'XTickLabel',scl.f_label);
        
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
    subplot(pp.sp_d(1)+1,pp.sp_d(2),subplot_dummy)
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
        subplot(pp.sp_d(1)+1,pp.sp_d(2),splot); caxis([c_diff(1) c_diff(2)]);
    else
        subplot(pp.sp_d(1)+1,pp.sp_d(2),splot); caxis([c(1) c(2)]);
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
            topoplot(ispc_data,chan_locs,'maplimits',[ispc_diff_limits(1) ispc_diff_limits(2)],'electrodes',pp.topo_elecs,'colormap',pp.cmap);
        else
            ispc_data=squeeze(mean(mean(abs(sum(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,cond,s_inds_g(:,group)),5))./ ...
                sum(abs(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,cond,s_inds_g(:,group))),5),3),1));
            %ispc_data=squeeze(mean(mean(abs(sum(exp(1i*angle(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,cond,s_inds_g(:,9)))),5))./ ...
            %    sum(abs(exp(1i*angle(wave_evkdata(t_start:t_end,pp.chosen_topochan,f_indiv,cond,s_inds_g(:,9))))),5),3),1));
            topoplot(ispc_data,chan_locs,'maplimits',[ispc_topo_scale(1) ispc_topo_scale(2)],'electrodes',pp.topo_elecs,'colormap',pp.cmap);
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
        %    linspace(phase_axes(1),phase_axes(2),100)*p_phase(1,cond,win)+p_phase(2,cond,win),scl.g_color{group});
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
            linspace(phase_axes(1),phase_axes(2),100)*p_phase(2,cond,win)+p_phase(1,cond,win),scl.g_color{group});
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