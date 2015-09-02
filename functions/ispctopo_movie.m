%% animate phase over time as a wind-rose

fs=256;
spacefactor=6;
outfile='MGT-Phase-Film.mat';
film_dum=0;
tf_scheme='cubicl'; %'linlhot' is also good
cmap=colormap(pmkmp(256,tf_scheme));


t_start_ms=[0 100 200 300 500 700];
t_end_ms=[100 200 300 500 700 900];
maxwins=length(t_start_ms);

f_start_hz=[5];
f_end_hz=[5];

conds_diff={[1 2],[3 4]};
chosen_groups=[9];

ispc_topo_scale=[0.06 0.80];
ispc_diff_limits=[-.12 .18];

%
tf_scheme='cubicl'; %'linlhot' is also good
cmap=colormap(pmkmp(256,tf_scheme));

v=zeros(maxconds+1,length(f_start_hz),length(t_start_ms),2);
figure;
for freq_range=1:length(f_start_hz)
    [c,f_start]=min(abs(freqs-f_start_hz(freq_range)));
    [c,f_end]=min(abs(freqs-f_end_hz(freq_range)));
    for cond=1:maxconds+1
        subplot(1,maxconds+1,subplot_dummy);
        for group=chosen_groups
        if cond==maxconds+1
            ispc_data=squeeze(mean(mean(mean(abs(sum(wavelet_evk_all(t_start:t_end,:,f_end:f_start,conds_diff{1},s_inds_g(:,9)),5))./ ...
                sum(abs(wavelet_evk_all(t_start:t_end,:,f_end:f_start,conds_diff{1},s_inds_g(:,9))),5),4),3),1)) - ...
                squeeze(mean(mean(mean(abs(sum(wavelet_evk_all(t_start:t_end,:,f_end:f_start,conds_diff{2},s_inds_g(:,9)),5))./ ...
                sum(abs(wavelet_evk_all(t_start:t_end,:,f_end:f_start,conds_diff{2},s_inds_g(:,9))),5),4),3),1));
            topoplot(ispc_data,chan_locs,'maplimits',[ispc_diff_limits(1) ispc_diff_limits(2)],'electrodes','off','colormap',cmap);
        else
            ispc_data=squeeze(mean(mean(abs(sum(wavelet_evk_all(t_start:t_end,:,f_end:f_start,cond,s_inds_g(:,9)),5))./ ...
                sum(abs(wavelet_evk_all(t_start:t_end,:,f_end:f_start,cond,s_inds_g(:,9))),5),3),1));
            topoplot(ispc_data,chan_locs,'maplimits',[ispc_topo_scale(1) ispc_topo_scale(2)],'electrodes','off','colormap',cmap);
        end
        v(cond,freq_range,win,1) = min(ispc_data); v(cond,freq_range,win,2) = max(ispc_data);
        title(sprintf('%s, %d - %d ms, %d - %d Hz',cond_label{cond},t_start_ms(win),t_end_ms(win),f_start_hz(freq_range),f_end_hz(freq_range)))
        end
    cbar
    end
end
tightfig;
cmin=min(min(min(v(1:end-1,:,:,1)))); cmax=max(max(max(v(1:end-1,:,:,2))));
cmin_diff=min(min(v(end,:,:,1))); cmax_diff=max(max(v(end,:,:,2)));


for chan=chosen_chans
for freq_range=1:length(f_start_hz);
%convert freqs to pts
[c,f_start]=min(abs(freqs-f_start_hz(freq_range)));
[c,f_end]=min(abs(freqs-f_end_hz(freq_range)));
for cond=1
    for timeframe=1:spacefactor:maxtimepts
        film_dum=film_dum+1;
        %calculate
        phase_data=rad2deg(squeeze(mean(mean(angle(wavelet_evk_all(timeframe, ...
            chan,f_end:f_start,cond,s_inds)),1),3)));
        intensity_data=squeeze(mean(mean(itc_results_all(timeframe, ...
            chan,f_end:f_start,cond,s_inds),1),3));
        %plot
        wind_title=sprintf('%d ms, %d Hz, %s',round(time_ms(timeframe)),...
            f_start_hz(freq_range),chan_label{chan});
        wind_rose(phase_data,intensity_data,'n',18,'di',[0.2:0.1:0.8],'ci',3,...
            'labtitle',wind_title,'lablegend','coherence','cmap',cmap,'quad',3);
        set(gcf,'position',[100 100 500 500]);
        h=gcf;
        film(film_dum)=getframe(h);
        close(h);
    end
end
end
end

save(outfile,'film')

%%

[a,w,p]=size(film(1).cdata);
h2=figure;
set(h2,'position',[100 100 w a]);
axis off
movie(h2,film,3,3)
%mov=immovie(film);

