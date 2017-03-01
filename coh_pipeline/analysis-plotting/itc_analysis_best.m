%% plot ERPS with groups superimposed

%warning off MATLAB:linkaxes:RequireDataAxes

%pp.figdum=pp.figdum_init;

% [20 22 24 26 28 49 51 59] %
% cyan = G1, magenta = G2
lims = [-8 19];
%lims = [-.25 .25];

h_line=zeros(length(pp.plotn_cond),length(pp.chosen_g));
for chan=[7 16 25] %pp.chosen_chan(pp.plotn_chan)
    pp.figdum=pp.figdum+1;
    figure;
    %figure(pp.figdum)
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


%% plot ERPs / topographies for individual subjects with mike cohen's ERPviewer

suspects = [5 6 20 23 29 33 39 56 60];

%subject = 1;

for s=52 %1:60 %suspects
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


%% single-subject significant ITC check figure
% butterfly plot for conditions and difference with major electrodes bolded
% GFP for conditions
% topographies at key time-point (~250 ms) for conditions and differences
% ERSP / ITC at key channel (Fz)
% ERCOH for key relation (intra-frontal)

s_inds=find(s_inds_g(:,3))';
%rejmat=false(size(s_inds));

%key channel
ck_chan=7;

%key regional hypothesis
ck_region=1;
plot_hypinds=find(opt.pair_inds==ck_region);

%baseline
[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

cmap=makecmap([-1 1],0);

% p-values for masking critical data
alpha = 0.05;

figure;
for s=s_inds
clear yn
clf
subplot_dummy=0;
for cond=pp.plotn_cond
    
    if cond==imp.maxconds+1
        itc_plot_data = bsxfun( @minus, squeeze(wave_evknormdata(:,ck_chan,:,pp.cond_diff{1},s)), ...
            squeeze(wave_evknormdata(:,ck_chan,:,pp.cond_diff{2},s)) );
%         ercoh_plot_data = bsxfun( @minus, meanx(cohdata(:,:,pp.cond_diff{1},plot_hypinds,s),[1 2]), ...
%             meanx(cohdata(:,:,pp.cond_diff{2},plot_hypinds,s),[1 2]) );
    else
        mask_val = sqrt(-log(alpha)/n_trials_all(cond, s));
        itc_plot_data = squeeze(wave_evknormdata(:,ck_chan,:,cond,s));
        itc_plot_data(itc_plot_data < mask_val) = 0;
%         ercoh_plot_data = bsxfun( @minus, meanx(cohdata(:,:,cond,plot_hypinds,s),[1 2]), ...
%             meanx(cohdata(t_start_b:t_end_b,:,:,plot_hypinds,s),2) ); %subtracting common baseline
    end
    
    %itc
    subplot_dummy=subplot_dummy+1; subplot(1,length(pp.plotn_cond),subplot_dummy);
    imagesc(fliplr(itc_plot_data)');
    set(gca, 'YDir', 'normal');
    % [~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour); set(h,'EdgeColor','None'); 
    axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms);
    set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label);
    grid on; set(gca,'Layer','Top'); hold on;
    plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
    colorbar;
    
    %ercoh
%     subplot_dummy=subplot_dummy+1; subplot(length(pp.plotn_cond),2,subplot_dummy);
%     [~,h]=contourf(fliplr(ercoh_plot_data)',pp.n_contour);
%     set(h,'EdgeColor','None'); axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
%     set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms);
%     set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label);
%     grid on; set(gca,'Layer','Top'); hold on;
%     plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
%     colorbar;
    
end

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

%erp_topo_scale=[-8 11];
erp_topo_scale=[-.2 .2];
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
for chan=[25] %pp.chosen_chan(pp.plotn_chan)
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


%% image PURE ITC in time-freq at a chosen channel, as increase from baseline
% add "with ERP superimposed"

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

p_alpha=0.05;
n_perms=2000;
do_statmask=true;

pp.figdum=pp.figdum_init;
%pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
clear overtitle
for group=pp.chosen_g(pp.plotn_g)
for chan=7 % 7 24 26 %pp.chosen_chan(pp.plotn_chan)
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
        itc_diff_data{cond}=bsxfun(@minus,squeeze(wave_evknormdata(:,chan,:,cond,s_inds_g(:,group))), ...
            itc_plot_base);
        %itc_diff_data{cond}=squeeze(wave_evknormdata(:,chan,:,cond,s_inds_g(:,group)));
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

%% image PURE ITC in time-freq at a chosen channel, as increase from baseline
% includes ability to mask condition data based on the critical ITPC value

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

p_alpha=0.05;
n_perms=2000;
do_statmask=true;

itc_crit = sqrt( bsxfun(@rdivide, -log(p_alpha), mean(n_trials_all,2) )) ;
custom_chanset = [9 7 8 27 25 26 54 58 55];

pp.figdum=pp.figdum_init;
%pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
clear overtitle
% for subject=1
for group=1 %pp.chosen_g(pp.plotn_g)
for chan=custom_chanset % 7 24 26 %pp.chosen_chan(pp.plotn_chan)
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
        itc_diff_data{cond}=bsxfun(@minus,squeeze(wave_evknormdata(:,chan,:,cond,s_inds_g(:,group))), ...
            itc_plot_base);
        %itc_diff_data{cond}=squeeze(mean(wave_evknormdata(:,chan,:,cond,:), 5));
        itc_plot_data = squeeze(mean( itc_diff_data{cond} , 3));
        itc_plot_data(itc_plot_data < itc_crit(cond)) = 0;
        imagesc(fliplr(itc_plot_data)');
        set(gca, 'YDir', 'normal');
        %[~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
        %set(h,'EdgeColor','None');
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
% end
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


%% image PURE ITC in time-freq at a chosen channel, as increase from baseline
% includes ability to mask condition data based on the critical ITPC value
% method 2: uses a comparison against a pre-stimulus baseline value

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
do_statmask=true;

itc_crit = sqrt( bsxfun(@rdivide, -log(p_alpha), mean(n_trials_all,2) )) ;
custom_chanset = [9 7 8 27 25 26 54 58 55];

pp.figdum=pp.figdum_init;
%pp.figdum_init=pp.figdum;
v=zeros(length(pp.chosen_g),length(pp.chosen_chan(pp.plotn_chan)),length(pp.plotn_cond),2);
clear overtitle
% for subject=1
for group=1 %pp.chosen_g(pp.plotn_g)
for chan=7 %custom_chanset % 7 24 26 %pp.chosen_chan(pp.plotn_chan)
pp.figdum=pp.figdum+1;
figure(pp.figdum); subplot_dummy=0;
overtitle{pp.figdum}=sprintf('%s / %s',scl.chan_label{chan},scl.g_label{group});
%itc_plot_base=permute(meanx(wave_evknormdata(t_start_b:t_end_b,chan,:,:,s_inds_g(:,group)),[3 5]),[3 1 2]);
%itc_plot_base=zeros(size(permute(meanx(wave_evknormdata(t_start_b:t_end_b,chan,:,:,s_inds_g(:,group)),[3 5]),[3 1 2])));
for cond=pp.plotn_cond
    subplot_dummy=subplot_dummy+1;
    subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
    if cond==imp.maxconds+1
        itc_plot_data= squeeze(mean( itc_diff_data{2} , 3));
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
        if cond==1 % baseline tiled
            itc_diff_data{cond}=repmat(permute(squeeze(mean(mean(...
                wave_evknormdata(t_start_b:t_end_b,chan,:,:,:), 1), 4)), ...
                [3 1 2]), ...
                [size(wave_evknormdata, 1) 1 1]);
        elseif cond==2 % actual condition mean data
            itc_diff_data{cond}=squeeze(mean(...
                wave_evknormdata(:,chan,:,:,:), 4));
        end
        itc_plot_data = squeeze(mean( itc_diff_data{cond} , 3));
        itc_plot_data(itc_plot_data < itc_crit(cond)) = 0;
        imagesc(fliplr(itc_plot_data)');
        set(gca, 'YDir', 'normal');
        %[~,h]=contourf(fliplr(itc_plot_data)',pp.n_contour);
        %set(h,'EdgeColor','None');
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
% end
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
for freq_range=1:2
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

%% image coherence in time-freq among a relational hypothesis' pairs (increase from common baseline)

sp_rowlabel={''};
sp_columnlabel=scl.cond_label;
x_plotlabel='Time (ms)';
y_plotlabel='Frequency (Hz)';

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

alpha=0.05;
n_perms=2000;

%pp.figdum_init=8;
pp.figdum_init=pp.figdum;
%pp.figdum=pp.figdum_init;
v=zeros(length(pp.chosen_g),6,length(pp.plotn_cond),2); clear overtitle;
for group=pp.chosen_g(pp.plotn_g)
for hyp=4:6
    pp.figdum=pp.figdum+1;
    figure(pp.figdum); subplot_dummy=0;
    if hyp==max(opt.pair_inds)+1
        overtitle{pp.figdum}='';
        plot_hypinds=1:imp.maxpairs;
    else
        overtitle{pp.figdum}=sprintf('%s / %s',opt.pair_indlbls{hyp},scl.g_label{group});
        plot_hypinds=find(opt.pair_inds==hyp);
    end
    coh_plot_base=permute(meanx(cohdata(t_start_b:t_end_b,:,:,plot_hypinds,s_inds_g(:,group)),[2 5]),[3 1 2]);
    for cond=pp.plotn_cond
        subplot_dummy=subplot_dummy+1;
        subplot(pp.sp_d(1),pp.sp_d(2),subplot_dummy)
        if cond==imp.maxconds+1
            coh_plot_data = meanx( bsxfun(@minus, ...
                cohdata(:,:,pp.cond_diff{1},plot_hypinds,s_inds_g(:,group)), ...
                cohdata(:,:,pp.cond_diff{2},plot_hypinds,s_inds_g(:,group))), [1 2]);
            [stats, df, pvals, surrog] = statcond( coh_diff_data, 'paired','on','method','perm', ...
                'naccu', n_perms, 'alpha', alpha, 'structoutput', 'on');
            %[p_fdr, p_masked] = fdr( stats.pval, alpha, 'nonParametric');
            [p_masked, p_fdr, adj_ci_cvrg, adj_p] = fdr_bh(pvals, alpha, 'pdep');
            v(group, hyp, cond,:) = [minx(coh_plot_data) maxx(coh_plot_data)];
            %coh_plot_data(~pvals) = 0;
            coh_plot_data(~p_masked) = 0;
            imagesc(fliplr(coh_plot_data)');
            set(gca,'YDir','normal');
            %[~,h]=contourf(fliplr(coh_plot_data)',pp.n_contour);
            %set(h,'EdgeColor','None');
        else
            coh_diff_data{cond}=bsxfun(@minus,meanx(cohdata(:,:,cond,plot_hypinds,s_inds_g(:,group)),[1 2 5]), ...
                coh_plot_base);
            %coh_diff_data{cond}=meanx(cohdata(:,:,cond,plot_hypinds,s_inds_g(:,group)),[1 2 5]);
            coh_plot_data = squeeze(mean( coh_diff_data{cond} , 3));
            [~,h]=contourf(fliplr(coh_plot_data)',pp.n_contour);
            set(h,'EdgeColor','None');
            v(group, hyp, cond,:) = caxis;
        end
        axis([scl.t_start scl.t_end 1 imp.maxfreqs]);
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms); %xlabel('Time (ms)');
        set(gca,'YTick',scl.f_ytick,'YTickLabel',scl.f_label); %ylabel('Frequency (Hz)');
        grid on; set(gca,'Layer','Top');
        hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(1,imp.maxfreqs,imp.maxfreqs),'k--'); hold off;
    end
    adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,pp.sp_d);    
end
end
v(v==0)=NaN;
c(1)=nanmin(nanmin(nanmin(v(:,:,1:end-1,1)))); c(2)=nanmax(nanmax(nanmax(v(:,:,1:end-1,2))));
c_diff(1)=nanmin(nanmin(nanmin(v(:,:,end,1)))); c_diff(2)=nanmax(nanmax(nanmax(v(:,:,end,2))));
cmap=makecmap(c,0,'purpleorange');
%cmap=parula(256);
cmap_diff=makecmap(c_diff,0,'redgreen');
for fig=pp.figdum_init+1:pp.figdum
    figure(fig);
    set(gcf,'position',[1300 120 800 300]);
    for splot=1:subplot_dummy
        temp_s=subplot(pp.sp_d(1),pp.sp_d(2),splot);
        if splot==subplot_dummy
            caxis([c_diff(1) c_diff(2)]);
            colormap(temp_s,cmap_diff);
        else
            caxis([c(1) c(2)]);
            colormap(temp_s,cmap);
        end
    end
    set_print_size(18,6);
    plottitle(overtitle{fig},2);
end
%make color bars separately
figure;
colorscale_plot(c, cmap, 0.25);
colorscale_plot(c_diff, cmap_diff, 0.75);
%
clear_plotassistvars

%% coherence in frequency band as topoplot with colored lines indicating strength, multiple time windows

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel=make_timelabels(pp.t_start_ms(1:end-1),pp.t_end_ms(1:end-1));
%sp_columnlabel{end}='';
x_plotlabel='Time Windows';
y_plotlabel='Conditions';

linewidth=3;
linescale=[1,256];
cbar_ticks=5;

line_limit=[-1 1];
line_limit_diff=[-.5 .5];

%line_limit=[-.02 .19];
%line_limit_diff=[-0.08 .08];

cmap=makecmap(line_limit,0,'purpleorange');
cmap_diff=makecmap(line_limit_diff,0,'redgreen');

n_perms = 2000;
p_alpha = .05;
plot_sig = true;

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

coh_linescale_mat=zeros(length(pp.plotn_cond),length(pp.chosen_g),length(pp.t_start_ms),length(pp.plotn_f),imp.maxpairs);
dummydata=ones(length(chan_locs),1)*0.3;
pp.figdum_init=pp.figdum;
%pp.figdum=pp.figdum_init;
%viz_pairs = 42:90;
viz_pairs = 1:41;
gdum2=0;
for group=pp.chosen_g(pp.plotn_g)
    gdum2=gdum2+1;
    fdum=0;
for freq_range=1:2 %pp.plotn_f
    fdum=fdum+1;
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    pp.figdum=pp.figdum+1;
    figure(pp.figdum); subplot_dummy=0; clear overtitle;
    overtitle{pp.figdum}=sprintf('Coherence Pairs, %s / %1.1f - %1.1f Hz',scl.g_label{group},pp.f_start_hz(freq_range),pp.f_end_hz(freq_range));
for cond=pp.plotn_cond
    for win=1:length(pp.t_start_ms)-2
        [~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
        [~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));
        subplot_dummy=subplot_dummy+1;
        subplot(length(pp.plotn_cond),length(pp.t_start_ms)-2,subplot_dummy)
        if win~=pp.maxwin
        topoplot(dummydata,chan_locs,'style','blank','maplimits',[0 1],'electrodes','off'); hold on;
        % do a significance test to determine which pairs should be plotted
        
        if cond==imp.maxconds+1 % if a difference, test difference
            statdata = cell(2,1);
            statdata{1}=meanx(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{1},:,s_inds_g(:,group)),[4 5]);
            statdata{2}=meanx(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{2},:,s_inds_g(:,group)),[4 5]);
            %[stats df pvals] = statcond( statdata, 'paired','on', 'structoutput', 'on');
            [stats df pvals] = statcond( statdata, 'paired','on', 'method', 'perm', 'naccu', n_perms, ...
                'alpha', p_alpha, 'structoutput', 'on');
            %[p_masked, p_fdr, adj_ci_cvrg, adj_p] = fdr_bh(pvals, p_alpha, 'pdep');
            %pairs2plot = find(p_masked)';
            pairs2plot = intersect(find(stats.mask),viz_pairs)';
            %pairs2plot = intersect(find(p_masked),viz_pairs)';
        else % if a condition, test increase from baseline
            statdata = cell(2,1);
            statdata{1}=meanx(cohdata(t_start:t_end,f_end:f_start,:,:,s_inds_g(:,group)),[4 5]);
            statdata{2}=meanx(cohdata(t_start_b:t_end_b,f_end:f_start,:,:,s_inds_g(:,group)),[4 5]);
            pairs2plot = viz_pairs;
        end
        
        for pair=pairs2plot
            %define [x1 x2], [y1 y2], and [z1 z2] of the arc
            x=[chan_locs(opt.coherence_pairs(pair,1)).topo_x chan_locs(opt.coherence_pairs(pair,2)).topo_x];
            y=[chan_locs(opt.coherence_pairs(pair,1)).topo_y chan_locs(opt.coherence_pairs(pair,2)).topo_y];
            %z=[2.1-randn*.05 2.1+randn*.05];
            %determine strength of coherence
            if cond==imp.maxconds+1
                paircoh_for_linecolor=mean(mean(mean(mean(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{1},pair,s_inds_g(:,group)),1),2),3),5) - ...
                    mean(mean(mean(mean(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{2},pair,s_inds_g(:,group)),1),2),3),5);
                coh_linescale_mat(cond,gdum2,win,fdum,pair)=paircoh_for_linecolor;
                paircoh_color=(paircoh_for_linecolor-line_limit_diff(1))/(line_limit_diff(2) - line_limit_diff(1))*(linescale(2)-1);
                paircoh_color = ceil(paircoh_color)+1;
                linecolor=cmap_diff(paircoh_color,:);
                linesize=norm2limits(paircoh_for_linecolor,line_limit_diff)*linewidth;
                linealpha=norm2limits(paircoh_for_linecolor,line_limit_diff);
            else
                paircoh_for_linecolor=mean(mean(mean(cohdata(t_start:t_end,f_end:f_start,cond,pair,s_inds_g(:,group)),1),2),5) -...
                    mean(mean(mean(mean(cohdata(t_start_b:t_end_b,f_end:f_start,:,pair,s_inds_g(:,group)),1),2),3),5); %subtract (common) baseline
                coh_linescale_mat(cond,gdum2,win,fdum,pair)=paircoh_for_linecolor;
                paircoh_color=(paircoh_for_linecolor-line_limit(1))/(line_limit(2) - line_limit(1))*(linescale(2)-1);
                paircoh_color = ceil(paircoh_color)+1;
                linecolor=cmap(paircoh_color,:);
                linesize=norm2limits(paircoh_for_linecolor,line_limit)*linewidth;
                linealpha=norm2limits(paircoh_for_linecolor,line_limit);
            end
            direction=1;
            plot_arc([x(1),y(1)], [x(2),y(2)], linecolor, linesize, direction);
            hold on;
        end
        else
            set(gca,'Visible','Off');
        end
        hold off;
    end
end
%adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle{fig},[length(pp.plotn_cond),length(pp.t_start_ms)-1], true);
set_print_size(16.5,16.5);
end
end
coh_linescale_mat(coh_linescale_mat==0)=NaN;
c=[nanmin(nanmin(nanmin(nanmin(nanmin(coh_linescale_mat(1:end-1,:,:,:,:)))))) nanmax(nanmax(nanmax(nanmax(nanmax(coh_linescale_mat(1:end-1,:,:,:,:))))))];
c_diff=[nanmin(nanmin(nanmin(nanmin(nanmin(coh_linescale_mat(end,:,:,:,:)))))) nanmax(nanmax(nanmax(nanmax(nanmax(coh_linescale_mat(end,:,:,:,:))))))];

figure;
colorscale_plot(line_limit, cmap, 0.25);
colorscale_plot(line_limit_diff, cmap_diff, 0.75);


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

%% baseline normalization options

%% channel subset options (can be regional)

%% scatter plots, error bar plots, boxplots, per-subject condition difference line plots

%% scatter baseline to event / condition difference in ITC / COH measures

%% correlation structure of frequency bands

%% things inside of coh_expfit

%% things inside of eec_analysis



