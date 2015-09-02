%% animate phase over time as a wind-rose

writerObj = VideoWriter('phasesum2.avi');
writerObj.FrameRate = 3;
open(writerObj);

outfile='PhaseSum-Film2.mat';

%%

f_start_hz=[4]; %3,5,8];
f_end_hz=[4]; %4,7,11];

maxconds=4;

conds_diff={[1 2],[3 4]};

tf_scheme='cubicl'; %'linlhot' is also good
cmap=colormap(pmkmp(256,tf_scheme));

chosen_chans=[56]; % 7 16 25];

itc_topo_scale=[0.1 0.5];
itc_diff_limits=[-.1 .05];
chosen_topochans=1:61;

di_bins=[0.1:0.1:0.7];

itc_lims=[0 0.4];

n_contours=9;

phase_axes=[-pi pi 0 1];

hist_axes=[0 1 0 15];
hist_nbins=10;

fs=256;
spacefactor=3;

%%

film_dum=0;
for chan=chosen_chans
for freq_range=1:length(f_start_hz);
%convert freqs to pts
[c,f_start]=min(abs(freqs-f_start_hz(freq_range)));
[c,f_end]=min(abs(freqs-f_end_hz(freq_range)));

%calculate time-invariant things
%itc_data=squeeze(mean(itc_results_all(:,chan,:,cond,s_inds_g(:,group)),5));

for timeframe=1:spacefactor:maxtimepts
    panel_phase;
    film_dum=film_dum+1;
    for cond=1
    gdum=0;
    %first populate the top left and right panels
    for group=chosen_groups
        gdum=gdum+1;
        %calculate things
        topo_data=squeeze(mean(mean(mean(itc_results_all(timeframe,...
            chosen_topochans,f_end:f_start,cond,s_inds_g(:,group)),1),3),5));
        phase_data=rad2deg(squeeze(mean(mean(angle(wavelet_evk_all(timeframe, ...
            chan,f_end:f_start,cond,s_inds_g(:,group))),1),3)));
        intensity_data=squeeze(mean(mean(itc_results_all(timeframe, ...
            chan,f_end:f_start,cond,s_inds_g(:,group)),1),3));
        itc_data=squeeze(mean(itc_results_all(:,chan,:,cond,s_inds_g(:,group)),5));
        erp_data=mean(mean_data_all(:,chan,cond,s_inds_g(:,group)),4);
        scatterdata_x=deunwrap(squeeze(unwrap(angle(wavelet_evk_all(timeframe, ...
            chan,f_start,cond,s_inds_g(:,group))))) - ...
            angle(mean(wavelet_evk_all(timeframe,chan,f_start,cond,s_inds_g(:,group)),5)));
        scatterdata_y=squeeze(itc_results_all(timeframe,chan,f_start,cond,s_inds_g(:,group)));
        itc_histdata=squeeze(mean(mean(mean(itc_results_all(timeframe,chan, ...
            f_end:f_start,cond,s_inds_g(:,group)),1),3),4));
        itc_histdata_se=std(squeeze(mean(mean(mean(itc_results_all(timeframe, ...
            chan,f_end:f_start,cond,s_inds_g(:,group)),1),3),4)));
        %[nels(:,group),cents(:,group)]=hist(itc_histdata,hist_nbins);
        [f(:,group),xi(:,group)]=ksdensity(itc_histdata);
        f(:,group)=f(:,group)/max(abs(f(:,group)))*hist_axes(4);
        %plot things
        %topo
        p(1,2*gdum-1).select();
        topoplot(topo_data,chan_locs,'maplimits',[itc_topo_scale(1) itc_topo_scale(2)],...
            'electrodes','off','colormap',cmap,'emarker2',{chan '.' g_color{group} 20 20});
        %windrose
        p(1,2*gdum).select();
        wrp=p(1,2*gdum).select();
        wind_rose(phase_data,intensity_data,'n',18,'ci',[10 20 30],...
            'lablegend','coherence','cmap',cmap,...
            'quad',3,'parent',wrp,'di',di_bins); %,'di',[0.1:0.1:0.7]);
        %ITC heatmap + ERP
        p(2,gdum).select();
        contourf(fliplr(itc_data)',n_contours); shading flat;
        caxis([itc_lims(1) itc_lims(2)]); colormap(cmap); hold on;
        plot(erp_data+8,'w'); hold on;
        vline(timeframe,'r-'); hold on; plot(timeframe,erp_data(timeframe)+8,'rx'); hold on;
        vline(time_zero,'k--'); hold on;
        hline(20-f_start,'w--'); hold off;
        axis([time_start_win time_end_win 1 maxfreqs]);
        set(gca,'XTick',time_xtick,'XTickLabel',time_xticks);
        set(gca,'YTick',freqs_ytick,'YTickLabel',freqs_label);
        grid on;
        % scatterplot
        p(3,1).select();
        scatter(scatterdata_x,scatterdata_y,g_color{group}); hold on;
        axis(phase_axes); grid on;
        %KDE
        p(3,2).select();
        plot(xi(:,group),f(:,group),g_color{group});
        hold on;
        plot(ones(1,100)*mean(itc_histdata),linspace(hist_axes(3),hist_axes(4),100),[g_color{group},'--']);
            hold on; plot(linspace(mean(itc_histdata)-itc_histdata_se,mean(itc_histdata)+itc_histdata_se,100), ...
            ones(1,100)*hist_axes(4)/20+gdum,g_color{group},'LineWidth',2); hold on;
        axis(hist_axes);
    end
    % histogram / KDEs
    %p(3,2).select();
    %bar(cents(:,chosen_groups),nels(:,chosen_groups),'stacked');
    %axis(hist_axes); hold on;
    %h=findobj(gca,'Type','patch');
    %set(h(1),'FaceColor',[0 1 0],'EdgeColor','w');
    %set(h(2),'FaceColor',[1 0 0],'EdgeColor','w');
    %plot(ones(1,100)*mean(itc_histdata),linspace(hist_axes(3),hist_axes(4),100),'k--'); hold on;
    %plot(linspace(mean(itc_histdata)-itc_histdata_se,mean(itc_histdata)+itc_histdata_se,100), ...
    %    ones(1,100)*hist_axes(4)/20,'k','LineWidth',2); hold off;
    end
    p(2,2).select(); set(gca,'YTickLabel',[]);
    set(gcf,'position',[0 120 1200 900]);
    h=gcf;
    film(film_dum)=getframe(h);
    writeVideo(writerObj,film(film_dum));
    close(h);
end
end
end

%%

save(outfile,'film')
close(writerObj);

%%

[a,w,p]=size(film(1).cdata);
h2=figure;
set(h2,'position',[0 120 w a]);
axis off
movie(h2,film,3,3)