%% animate phase over time as a wind-rose

writerObj = VideoWriter('windrose2.avi');
writerObj.FrameRate = 3;
open(writerObj);

outfile='MGT-Phase-Film.mat';

%%

f_start_hz=[5]; %3,5,8];
f_end_hz=[5]; %4,7,11];

maxconds=4;

conds_diff={[1 2],[3 4]};

tf_scheme='cubicl'; %'linlhot' is also good
cmap=colormap(pmkmp(256,tf_scheme));

chosen_chans=[25]; % 7 16 25];

fs=256;
spacefactor=6;

film_dum=0;
for chan=chosen_chans
for freq_range=1:length(f_start_hz);
figure; subplot_dummy=0;
%convert freqs to pts
[c,f_start]=min(abs(freqs-f_start_hz(freq_range)));
[c,f_end]=min(abs(freqs-f_end_hz(freq_range)));
for cond=1:maxconds+1
    for timeframe=1:spacefactor:maxtimepts
        film_dum=film_dum+1;
        %calculate
        if cond==maxconds+1
            phase_data=rad2deg(squeeze(...
                mean(mean(mean(angle(wavelet_evk_all(timeframe,chan,f_end:f_start,conds_diff{1},s_inds)),1),3),4) -...
                mean(mean(mean(angle(wavelet_evk_all(timeframe,chan,f_end:f_start,conds_diff{2},s_inds)),1),3),4)));
            intensity_data=squeeze(...
                mean(mean(mean(itc_results_all(timeframe,chan,f_end:f_start,conds_diff{1},s_inds),1),3),4) - ...
                mean(mean(mean(itc_results_all(timeframe,chan,f_end:f_start,conds_diff{2},s_inds),1),3),4));
        else
            phase_data=rad2deg(squeeze(mean(mean(angle(wavelet_evk_all(timeframe, ...
                chan,f_end:f_start,cond,s_inds)),1),3)));
            intensity_data=squeeze(mean(mean(itc_results_all(timeframe, ...
                chan,f_end:f_start,cond,s_inds),1),3));
        end
        %plot
        wind_title=sprintf('%d ms, %d Hz, %s',round(time_ms(timeframe)),...
            f_start_hz(freq_range),chan_label{chan});
        wind_rose(phase_data,intensity_data,'n',18,'di',[0.2:0.1:0.8],'ci',3,...
            'labtitle',wind_title,'lablegend','coherence','cmap',cmap,'quad',3);
        set(gcf,'position',[100 100 500 500]);
        h=gcf;
        film(film_dum)=getframe(h);
        writeVideo(writerObj,film(film_dum));
        close(h);
    end
end
end
end

save(outfile,'film')
close(writerObj);

%%

[a,w,p]=size(film(1).cdata);
h2=figure;
set(h2,'position',[100 100 w a]);
axis off
movie(h2,film,3,3)
%mov=immovie(film);


%%


