%% fix up some scaling / plotting stuff

%pp.plotn_cond=1:4;
%pp.sp_d = numSubplots(length(pp.plotn_cond));

[~, scl.t_start] = min(abs(scl.t_ms+500));
[~, scl.t_end] = min(abs(scl.t_ms-800));

%% plot ERPS with groups superimposed

%warning off MATLAB:linkaxes:RequireDataAxes

pp.figdum=pp.figdum_init;

% [20 22 24 26 28 49 51 59] %
% cyan = G1, magenta = G2
lims = [-4 6];

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
        axis([1 imp.erpmaxtimepts lims(1) lims(2)]);
        %axis tight
        vline(scl.t_zero_erp,'k--'); hold off;
        set(gca,'XTick',scl.t_xtick_erp,'XTickLabel',scl.t_xtick_ms)
        grid on
        title([scl.chan_label{chan},'/',scl.cond_label{cond}])
    end
tightfig;
%linkaxes(sp(1:end-1))
set_print_size(17,17/scl.phi);
set(gcf,'position',[120 120 987 300]);
end
clear_plotassistvars

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
for chan=25 %pp.chosen_chan(pp.plotn_chan)
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