
do_film=true;

if do_film

writerObj = VideoWriter('ern_cohlinetopowithERP_theta_nki_0001.avi');
writerObj.FrameRate = 4;
open(writerObj);

outfile='ern_cohlinetopowithERP_theta_nki_0001.mat';

end

%%

t_start_ms=-100;
t_end_ms=700;

spacefactor=1;

%f_start_hz=[2.9, 4,   8]; %3,5,8];
%f_end_hz=  [3, 5.3, 12.8]; %4,7,11];

f_start_hz=[2,   3.5, 4.2, 5.5,  8,    14,   27, 3, 5];     %lower edges
f_end_hz=[  3.2, 4.2, 5,   7.3,  13,   24,   32, 5, 7];  %upper edges

linescale=[1,256];
%line_limit=[-.09 .09];
%line_limit_diff=[-.06 .06];

line_limit=[-.04 .18];
line_limit_diff=[-.06 .13];

cmap=makecmap(line_limit);
cmap_diff=makecmap(line_limit_diff);

alpha=.001;
plot_sig=true;

erpchan=7;
erpgroup=9;
%

sp_columnlabel={scl.cond_label{pp.plotn_cond}};

x_plotlabel=' ';
y_plotlabel=' ';
subplot_dims=[3,3];
subplot_key=[1 4 7 2 5 8 3 6 9];

cb_pos=0.25;
cb_width=0.02;
n_colors=256;
n_ticks=5;

warning off MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0

%%

film_dum=0;

[~,t_start]=min(abs(scl.t_ms-t_start_ms));
[~,t_end]=min(abs(scl.t_ms-t_end_ms));

dummydata=ones(length(chan_locs),1)*0.3;

for timeframe=t_start:spacefactor:t_end
    figure;
    set(gcf,'position',[0 120 1280 800]);
    film_dum=film_dum+1;
    overtitle=[num2str(round(scl.t_ms(timeframe))),' ms'];
    subplot_dummy=0;
    
    for freq_range=3
        
        sp_rowlabel={['CTL ',num2str(f_start_hz(freq_range)),'-',num2str(f_end_hz(freq_range)),' Hz'], ...
            ['ALC ',num2str(f_start_hz(freq_range)),'-',num2str(f_end_hz(freq_range)),' Hz'],'ERP'};
        %convert freqs to pts
        [~,f_start]=min(abs(scl.freqs-f_start_hz(freq_range)));
        [~,f_end]=min(abs(scl.freqs-f_end_hz(freq_range)));
        
    for cond=pp.plotn_cond
    for group=pp.chosen_g(pp.plotn_g)
    
    subplot_dummy=subplot_dummy+1;
    subplot(subplot_dims(1),subplot_dims(2),subplot_key(subplot_dummy));
    
    %
    topoplot(dummydata,chan_locs,'style','blank','maplimits',[0 1],'electrodes','off'); hold on;
    for pair=1:imp.maxpairs
        %determine significance before further plotting
        if plot_sig
        gdum=0;
        for statsgroup=pp.chosen_g(pp.plotn_g)
            gdum=gdum+1;
            ranova_data{gdum}=squeeze(mean(mean(cohdata(timeframe,f_end:f_start,:,pair,s_inds_g(:,statsgroup)),1),2))';
        end
        [p,table]=anova_rm(ranova_data,'off');
        %if p(1) > alpha
        %if p(2)>alpha
        %if p(4) > alpha
        if ~any(p([1:2,4])<alpha)
            continue
        end
        end
        %define [x1 x2], [y1 y2], and [z1 z2] of the arc
        x=[chan_locs(opt.coherence_pairs(pair,1)).topo_x chan_locs(opt.coherence_pairs(pair,2)).topo_x];
        y=[chan_locs(opt.coherence_pairs(pair,1)).topo_y chan_locs(opt.coherence_pairs(pair,2)).topo_y];
        %determine strength of coherence
        if cond==imp.maxconds+1
            paircoh_for_linecolor=mean(mean(mean(mean(cohdata(timeframe,f_end:f_start,pp.cond_diff{1},pair,s_inds_g(:,group)),1),2),3),5) - ...
                mean(mean(mean(mean(cohdata(timeframe,f_end:f_start,pp.cond_diff{2},pair,s_inds_g(:,group)),1),2),3),5);
            paircoh_color=(paircoh_for_linecolor-line_limit_diff(1))/(line_limit_diff(2) - line_limit_diff(1))*(linescale(2)-1);
            paircoh_color = ceil(paircoh_color)+1;
            linecolor=cmap_diff(paircoh_color,:);
            linesize=norm2limits(paircoh_for_linecolor,line_limit_diff)*5;
            linealpha=norm2limits(paircoh_for_linecolor,line_limit_diff);
        else
            paircoh_for_linecolor=mean(mean(mean(cohdata(timeframe,f_end:f_start,cond,pair,s_inds_g(:,group)),1),2),5) -...
                mean(mean(mean(mean(cohdata(1:scl.t_zero,f_end:f_start,:,pair,s_inds_g(:,group)),1),2),3),5); %subtract baseline
            paircoh_color=(paircoh_for_linecolor-line_limit(1))/(line_limit(2) - line_limit(1))*(linescale(2)-1);
            paircoh_color = ceil(paircoh_color)+1;
            linecolor=cmap(paircoh_color,:);
            linesize=norm2limits(paircoh_for_linecolor,line_limit)*5;
            linealpha=norm2limits(paircoh_for_linecolor,line_limit);
        end
        direction=1;
        %plot the arc
        patchline_arc([x(1),y(1)], [x(2),y(2)], linecolor, linesize, linealpha, direction);
        hold on;
    end
    hold off
    if cond==imp.maxconds+1
        colormap(cmap_diff)
        cb_ax(cb_pos,cb_width,line_limit_diff,n_colors,n_ticks,'vert');
    else
        colormap(cmap)
        cb_ax(cb_pos,cb_width,line_limit,n_colors,n_ticks,'vert');
    end
    freezeColors;
    
    end
    
    %plot an ERP on the side
    if true
    subplot_dummy=subplot_dummy+1;
    subplot(subplot_dims(1),subplot_dims(2),subplot_key(subplot_dummy))
    if cond==imp.maxconds+1
        erp_plot_data=meanx(erpdata(:,erpchan,pp.cond_diff{1},s_inds_g(:,erpgroup)),1) -...
            meanx(erpdata(:,erpchan,pp.cond_diff{2},s_inds_g(:,erpgroup)),1);
        erp_plot_data_std = std(mean(erpdata(:,erpchan,pp.cond_diff{1},s_inds_g(:,erpgroup)),3) - ...
            mean(erpdata(:,erpchan,pp.cond_diff{2},s_inds_g(:,erpgroup)),3),0,4)/sqrt(sum(s_inds_g(:,erpgroup)));
    else
        erp_plot_data=meanx(erpdata(:,erpchan,cond,s_inds_g(:,erpgroup)),1);
        erp_plot_data_std=std(erpdata(:,erpchan,cond,s_inds_g(:,erpgroup)),0,4)/sqrt(sum(s_inds_g(:,erpgroup)));
    end
    if false
        shadedErrorBar(1:imp.maxtimepts,erp_plot_data,erp_plot_data_std); hold on;
        plot(ones(100,1)*timeframe,linspace(-15,20,100),'r'); hold on;
        axis([scl.t_start scl.t_end -.15 .20]);
        vline(scl.t_zero,'k--'); hold off;
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms);
        xlabel('Time (ms)');
        %ylabel('Amplitude (uV/cm^2)');
        grid on;
    else
        plot(erp_plot_data); hold on;
        plot(ones(100,1)*timeframe,linspace(-.15,.20,100),'r'); hold on;
        axis([scl.t_start scl.t_end -.15 .20]);
        vline(scl.t_zero,'k--'); hold off;
        set(gca,'Visible','Off'); %it's possible that "axis off" does this
    end
    end
    end
    makescale([140 .10],.05,round(.2*opt.timerate),'200');
    end
    adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,subplot_dims);
    %tightfig;
    %set(gcf,'position',[0 120 1280 800]);
    h=gcf;
    if do_film
        
        film(film_dum)=getframe(h);
        writeVideo(writerObj,film(film_dum));
        
    end
    close(h);
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