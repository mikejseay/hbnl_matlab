do_film=true;

if do_film

writerObj = VideoWriter('err_erptopo_947.avi');
writerObj.FrameRate = 3;
open(writerObj);

outfile='err_erptopo_947.mat';

end

%%

t_start_ms=-100;
t_end_ms=700;

erp_topo_scale=[-16 22];
erp_diff_limits=[-9 6];

erpchan=7;

spacefactor=1;

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

cmap=makecmap(erp_topo_scale);
cmap_diff=makecmap(erp_diff_limits);

sp_rowlabel={'Topography','ERP at FZ'}; %{scl.cond_label{pp.plotn_cond}};
sp_columnlabel={scl.cond_label{pp.plotn_cond}}; %{scl.g_label{pp.chosen_g}};
x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[2, 3];
subplot_key=[1 4 2 5 3 6];

%%
%film=zeros(length(1:spacefactor:imp.maxtimepts),1);
film_dum=0;
%for freq_range=1:length(f_start_hz);

%convert freqs to pts
%[~,f_start]=min(abs(scl.freqs-f_start_hz(freq_range)));
%[~,f_end]=min(abs(scl.freqs-f_end_hz(freq_range)));

[~,t_start]=min(abs(scl.t_ms-t_start_ms));
[~,t_end]=min(abs(scl.t_ms-t_end_ms));

%calculate time-invariant things
%itc_data=squeeze(mean(itc_results_all(:,chan,:,cond,s_inds_g(:,group)),5));

for timeframe=t_start:spacefactor:t_end
    figure;
    film_dum=film_dum+1;
    overtitle=[num2str(round(scl.t_ms(timeframe))),' ms'];
    subplot_dummy=0;
    for cond=pp.plotn_cond
    for group=pp.chosen_g
    subplot_dummy=subplot_dummy+1;
    subplot(subplot_dims(1),subplot_dims(2),subplot_key(subplot_dummy));
    %calculate things
    if cond==imp.maxconds+1
    topo_data=meanx(erpdata(timeframe,:,pp.cond_diff{1},s_inds_g(:,group)),2) -...
        meanx(erpdata(timeframe:timeframe+spacefactor,:,pp.cond_diff{2},s_inds_g(:,group)),2);
    topoplot(topo_data',chan_locs,'maplimits',erp_diff_limits,...
        'electrodes','off','colormap',cmap_diff,'style','fill','numcontour',pp.n_contour);
    freezeColors;
    else
    topo_data=meanx(erpdata(timeframe:timeframe+spacefactor,:,cond,s_inds_g(:,group)),2);
    topoplot(topo_data',chan_locs,'maplimits',erp_topo_scale,...
        'electrodes','off','colormap',cmap,'style','fill','numcontour',pp.n_contour);
    freezeColors;
    end
    shading flat
    end
    %colorbar
    if false
        colorbar;
        cbfreeze;
    else
        cb_pos = get(gca,'Position') + [0.24 0 0 0];
        cb_pos(3) = 0.03; %set width
        if cond==imp.maxconds+1
        cb_lim = erp_diff_limits;
        else
        cb_lim = erp_topo_scale;
        end
        colorscale([1 256], cb_lim, range(cb_lim)/5, 'vert', ...
            'Position',cb_pos);
        freezeColors;
    end
        
    %erp part
    subplot_dummy=subplot_dummy+1;
    subplot(subplot_dims(1),subplot_dims(2),subplot_key(subplot_dummy));
    if cond==imp.maxconds+1
        erp_plot_data=meanx(erpdata(:,erpchan,pp.cond_diff{1},s_inds_g(:,group)),1) -...
            meanx(erpdata(:,erpchan,pp.cond_diff{2},s_inds_g(:,group)),1);
        erp_plot_data_std = std(mean(erpdata(:,erpchan,pp.cond_diff{1},s_inds_g(:,group)),3) - ...
            mean(erpdata(:,erpchan,pp.cond_diff{2},s_inds_g(:,group)),3),0,4)/sqrt(sum(s_inds_g(:,group)));
    else
        erp_plot_data=meanx(erpdata(:,erpchan,cond,s_inds_g(:,group)),1);
        erp_plot_data_std=std(erpdata(:,erpchan,cond,s_inds_g(:,group)),0,4)/sqrt(sum(s_inds_g(:,group)));
    end
    plot(erp_plot_data); hold on;
    plot(ones(100,1)*timeframe,linspace(-15,20,100),'r'); hold on;
    axis([scl.t_start scl.t_end -15 20]);
    vline(scl.t_zero,'k--'); hold off;
    %set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms);
    %xlabel('Time (ms)');
    %ylabel('Amplitude (uV/cm^2)');
    %grid on;
    %box off;
    set(gca,'Visible','Off');
    end
    makescale([140 10],5,round(.2*param_struct.timerate),'200');
    adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
    %subplot(3,2,subplot_dummy); colorbar; cbfreeze;
    if do_film
        set(gcf,'position',[0 120 1280 800]);
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