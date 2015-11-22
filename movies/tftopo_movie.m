%% animate phase over time as a wind-rose

do_film=true;

if do_film

writerObj = VideoWriter('ern_tftopo_3freqs_947.avi');
writerObj.FrameRate = 3;
open(writerObj);

outfile='ern_tftopo_3freqs_947.mat';

end

%%

t_start_ms=-100;
t_end_ms=700;

spacefactor=1;

f_start_hz=[2.9, 4,   8]; %3,5,8];
f_end_hz=  [3,   5.3, 12.8]; %4,7,11];


%ero_topo_scale=[-1 15;-1 25;-25 8];
%ero_diff_limits=[-3 3;-0.8 10;-0.9 2.5];
ero_topo_scale=[-25 25;-25 25;-25 25];
ero_diff_limits=[-3 10;-3 10;-3 10];

itc_topo_scale=[0.2 0.5; 0.2 0.5;0.2 0.5];
itc_diff_limits=[-0.06 0.08; -.06 .08; -.06 0.08];
%itc_topo_scale=[0.2 0.45; 0.2 0.5;0.1 0.4];
%itc_diff_limits=[-0.06 0.03; -.03 .08; -.04 0.04];
%itc_topo_scale=[0.1 0.4; 0.1 0.5;0.1 0.4];
%itc_diff_limits=[-0.08 0.06; -.05 .12; -.06 0.1];

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel={'ERO-Delta','ITC-Delta','ERO-Theta','ITC-Theta','ERO-Alpha','ITC-Alpha'};
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[3,6];
%subplot_key=reshape(1:18,[6 3])';
subplot_key=reshape(1:18,[2 3 3]);

cb_pos=0.12;
cb_width=0.01;
n_colors=256;
n_ticks=5;

%%
%film=zeros(length(1:spacefactor:imp.maxtimepts),1);
film_dum=0;
%for freq_range=1:length(f_start_hz);

[~,t_start]=min(abs(scl.t_ms-t_start_ms));
[~,t_end]=min(abs(scl.t_ms-t_end_ms));

%calculate time-invariant things
%itc_data=squeeze(mean(itc_results_all(:,chan,:,cond,s_inds_g(:,group)),5));

for timeframe=t_start:spacefactor:t_end
    figure;
    set(gcf,'position',[80 120 1600 800]);
    film_dum=film_dum+1;
    overtitle=[num2str(round(scl.t_ms(timeframe))),' ms'];
    for freq_range=1:3
    %convert freqs to pts
    [~,f_start]=min(abs(scl.freqs-f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-f_end_hz(freq_range)));
    
    
    ero_cmap=makecmap(ero_topo_scale(freq_range,:));
    ero_cmap_diff=makecmap(ero_diff_limits(freq_range,:));
    itc_cmap=makecmap(itc_topo_scale(freq_range,:));
    itc_cmap_diff=makecmap(itc_diff_limits(freq_range,:));
    for cond=pp.plotn_cond
    for group=pp.chosen_g(1)
    subplot(subplot_dims(1),subplot_dims(2),subplot_key(1,freq_range,cond));
    %ero part
    if cond==imp.maxconds+1
    ero_topo_data=meanx(wave_totdata(timeframe,:,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),2) -...
        meanx(wave_totdata(timeframe:timeframe+spacefactor,:,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),2);
    topoplot(ero_topo_data',chan_locs,'maplimits',ero_diff_limits(freq_range,:),...
        'electrodes','off','colormap',ero_cmap_diff,'style','fill','numcontour',pp.n_contour);
    freezeColors;
    else
    ero_topo_data=meanx(wave_totdata(timeframe,:,f_end:f_start,cond,s_inds_g(:,group)),2) -...
        meanx(wave_totdata(1:scl.t_zero,:,f_end:f_start,cond,s_inds_g(:,group)),2);
    topoplot(ero_topo_data',chan_locs,'maplimits',ero_topo_scale(freq_range,:),...
        'electrodes','off','colormap',ero_cmap,'style','fill','numcontour',pp.n_contour);
    freezeColors;
    end
    shading flat
    %colorbar
    if cond==imp.maxconds+1
        cb_ax(cb_pos,cb_width,ero_diff_limits(freq_range,:),n_colors,n_ticks,'vert');
    else
        cb_ax(cb_pos,cb_width,ero_topo_scale(freq_range,:),n_colors,n_ticks,'vert');
    end
    freezeColors;
    %
    subplot(subplot_dims(1),subplot_dims(2),subplot_key(2,freq_range,cond));
    %itc part
    if cond==imp.maxconds+1
    itc_topo_data=meanx(itcdata(timeframe,:,f_end:f_start,pp.cond_diff{1},s_inds_g(:,group)),2) -...
        meanx(itcdata(timeframe:timeframe+spacefactor,:,f_end:f_start,pp.cond_diff{2},s_inds_g(:,group)),2);
    topoplot(itc_topo_data',chan_locs,'maplimits',itc_diff_limits(freq_range,:),...
        'electrodes','off','colormap',itc_cmap_diff,'style','fill','numcontour',pp.n_contour);
    freezeColors;
    else
    itc_topo_data=meanx(itcdata(timeframe,:,f_end:f_start,cond,s_inds_g(:,group)),2);
    topoplot(itc_topo_data',chan_locs,'maplimits',itc_topo_scale(freq_range,:),...
        'electrodes','off','colormap',itc_cmap,'style','fill','numcontour',pp.n_contour);
    freezeColors;
    end
    shading flat
    %colorbar
    if cond==imp.maxconds+1
        cb_ax(cb_pos,cb_width,itc_diff_limits(freq_range,:),n_colors,n_ticks,'vert');
    else
        cb_ax(cb_pos,cb_width,itc_topo_scale(freq_range,:),n_colors,n_ticks,'vert');
    end
    freezeColors;
    end
    end
    end
    adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
    %subplot(3,2,subplot_dummy); colorbar; cbfreeze;
    if do_film
        h=gcf;
        film(film_dum)=getframe(h);
        writeVideo(writerObj,film(film_dum));
        close(h);
    end
end
%end


%%

if do_film

save(outfile,'film')
close(writerObj);

%%

[a,w,p]=size(film(1).cdata);
h2=figure;
set(h2,'position',[80 120 w a]);
axis off
movie(h2,film,3,3)

end

clear_plotassistvars
clear_movieassistvars