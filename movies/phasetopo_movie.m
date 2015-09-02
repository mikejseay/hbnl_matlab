%% animate phase over time as a wind-rose

recording=false;
if recording 
writerObj = VideoWriter('phasetopo.avi');
writerObj.FrameRate = 3;
open(writerObj);

outfile='PhaseTopo-Film.mat';
end
%%

f_start_hz=[4]; %3,5,8];
f_end_hz=[4]; %4,7,11];

maxconds=4;

conds_diff={[1 2],[3 4]};

tf_scheme='cubicl'; %'linlhot' is also good
cmap=colormap(pmkmp(256,tf_scheme));

chosen_chans=[56]; % 7 16 25];

phase_topo_scale=[-pi pi];
phase_diff_limits=[-pi pi];

chosen_topochans=1:61;

fs=256;
spacefactor=3;

%%

film_dum=0;

%convert freqs to pts
[c,f_start]=min(abs(freqs-f_start_hz(freq_range)));
[c,f_end]=min(abs(freqs-f_end_hz(freq_range)));
for timeframe=1:spacefactor:maxtimepts
    %figure;
    film_dum=film_dum+1;
    subplot_dummy=0;
    for group=chosen_groups
    for cond=1:maxconds+1
        subplot_dummy=subplot_dummy+1;
        subplot(length(chosen_groups),maxconds+1,subplot_dummy)
        %calculate things
        if cond==maxconds+1
            topo_data=squeeze(angle(mean(mean(mean(mean(wavelet_evk_all(timeframe,:,f_end:f_start,conds_diff{1},s_inds_g(:,group)),1),3),4),5) -...
                mean(mean(mean(mean(wavelet_evk_all(t_end:t_end,:,f_end:f_start,conds_diff{2},s_inds_g(:,group)),1),3),4),5)));
            topoplot(contour_wrapedges(topo_data),chan_locs,'maplimits',[phase_diff_limits(1) phase_diff_limits(2)],'electrodes','off','colormap',cmap);
        else
            topo_data=squeeze(angle(mean(mean(mean(wavelet_evk_all(timeframe,:,f_end:f_start,cond,s_inds_g(:,group)),1),3),5)));
            topoplot(contour_wrapedges(topo_data),chan_locs,'maplimits',[phase_topo_scale(1) phase_topo_scale(2)],'electrodes','off','colormap',cmap);
        end
        %add group/condition labels
        if cond==1
            hold on;
            text(-1,0,g_label{group}); hold off;
        end
        if group==chosen_groups(1)
            hold on;
            text(-0.2,0.8,cond_label{cond}); hold off;
        end
    end
    end
    %timebar
    a=axes('Position', [0.1 0.5 0.8 0.1], 'Visible', 'on');
    set(a,'XLim',[1 maxtimepts],'XTick',time_xtick,'XTickLabel',time_xticks);
    xlabel('Time (ms)');
    %colorbar
    a=axes('Position', [0.88 0.15 0.12 0.7], 'Visible', 'off');
    set(a,'CLim',[-pi pi]);
    c=colorbar('YLim',[-pi pi],'YTick',[-pi,0,pi],'YTickLabel',{'0',[char(177),'pi/2'],[char(177),'pi']});
    %size and position the figure
    set(gcf,'position',[0 120 1000 600]);
    %take frames
    if recording
        h=gcf;
        film(film_dum)=getframe(h);
        writeVideo(writerObj,film(film_dum));
        close(h);
    end
end

%%

if recording
save(outfile,'film')
close(writerObj);

[a,w,p]=size(film(1).cdata);
h2=figure;
set(h2,'position',[0 120 w a]);
axis off
movie(h2,film,3,3)
end