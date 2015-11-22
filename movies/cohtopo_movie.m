do_film=true;

if do_film

writerObj = VideoWriter('ern_cohtopo_3freqs2.avi');
writerObj.FrameRate = 3;
open(writerObj);

outfile='ern_cohtopo_3freqs2.mat';

end

%%

t_start_ms=-100;
t_end_ms=700;

spacefactor=1;

f_start_hz=[2.9, 4,   8]; %3,5,8];
f_end_hz=  [3,   5.3, 12.8]; %4,7,11];

coh_topo_scale=[0 .12; 0 .15; -.08 .08];
coh_diff_limits=[-.06 .06; -.05 .05; -.04 .04];

%seed_pairs=[1:30;31:60;61:90;91:120];
%seed_label={'F4','F3','P3','P4'};
%seed_inds=[8,9,23,24]; %17 18
seed_pairs=[1:30;91:120];
seed_label={'F4','P4'};
seed_inds=[8,24]; %17 18

seed_chanlocs=load('/export/home/mike/matlab/origin/coords/31chans_ns.mat');
seed_chanlocs=getfield(seed_chanlocs,'chan_locs');

sp_rowlabel={scl.cond_label{pp.plotn_cond}};
sp_columnlabel={'F4-Delta','F4-Theta','F4-Alpha','P4-Delta','P4-Theta','P4-Alpha'};
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[3,6];
%subplot_key=reshape(1:18,[6 3])';
%subplot_key=reshape(1:18,[2 3 3]);

cb_pos=0.12;
cb_width=0.01;
n_colors=256;
n_ticks=5;

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

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
    subplot_dummy=0;
    
    for cond=pp.plotn_cond
    for group=pp.chosen_g(pp.plotn_g)
    for seed=1:2
    for freq_range=1:3
    subplot_dummy=subplot_dummy+1;
        
    %convert freqs to pts
    [~,f_start]=min(abs(scl.freqs-f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-f_end_hz(freq_range)));
    coh_cmap=makecmap(coh_topo_scale(freq_range,:));
    coh_cmap_diff=makecmap(coh_diff_limits(freq_range,:));
    
    subplot(subplot_dims(1),subplot_dims(2),subplot_dummy);
    %ero part
    if cond==imp.maxconds+1
    coh_topo_data=meanx(cohdata(timeframe,f_end:f_start,pp.cond_diff{1},seed_pairs(seed,:),...
        s_inds_g(:,group)),4) - ...
        meanx(cohdata(timeframe,f_end:f_start,pp.cond_diff{2},seed_pairs(seed,:),...
        s_inds_g(:,group)),4);
    coh_topo_data=insertrows(coh_topo_data,0,seed_inds(seed)-1);
    topoplot(coh_topo_data',seed_chanlocs,'maplimits',coh_diff_limits(freq_range,:),...
        'electrodes','off','colormap',coh_cmap_diff,'style','fill','numcontour',pp.n_contour);
    freezeColors;
    else
    coh_topo_data=meanx(cohdata(timeframe,f_end:f_start,cond,seed_pairs(seed,:),...
        s_inds_g(:,group)),4) -...
        meanx(cohdata(1:scl.t_zero,f_end:f_start,cond,seed_pairs(seed,:),...
        s_inds_g(:,group)),4);
    coh_topo_data=insertrows(coh_topo_data,0,seed_inds(seed)-1);
    topoplot(coh_topo_data',seed_chanlocs,'maplimits',coh_topo_scale(freq_range,:),...
        'electrodes','off','colormap',coh_cmap,'style','fill','numcontour',pp.n_contour);
    freezeColors;
    end
    shading flat
    %colorbar
    if cond==imp.maxconds+1
        cb_ax(cb_pos,cb_width,coh_diff_limits(freq_range,:),n_colors,n_ticks,'vert');
    else
        cb_ax(cb_pos,cb_width,coh_topo_scale(freq_range,:),n_colors,n_ticks,'vert');
    end
    freezeColors;
    end
    end
    end
    end
    adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
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