% mgt_coh_analysis.m - Analyses coherence results from the MGT task.

%% load in file list, parameter info, channel locations, and detect os

%pair_subset=[1:6 9 10 14]; %theta

load('/export/home/mike/matlab/coords/31chans_ns.mat')
load('/export/home/mike/matlab/mgt_coh/subjectid_table_ctl14.mat')
load('/export/home/mike/matlab/mgt_coh/ern_params_07_20_15.mat')

%% load in subject data

%declare number of subjects to look at
ns=100;

%size of data dimensions
%maxtimepts=round(param_struct.prestim_ms*param_struct.rate/1000)+param_struct.n_samps+1;
maxtimepts=param_struct.n_samps+1;
maxfreqs=length(param_struct.scale);
maxconds=length(param_struct.case_vec);
maxchans=length(param_struct.chan_vec);
maxpairs=size(param_struct.coherence_pairs,1);

%check if statistical information is present
coh_stats_present=param_struct.perms>0;

%pre-allocate vars
coh_results_all=zeros(maxtimepts,maxfreqs,maxconds,maxpairs,2);
coh_results_stats=zeros(maxtimepts,maxfreqs,maxconds,maxpairs,2);
mean_data_all=zeros(maxtimepts,maxchans,maxconds,2);
n_trials_all=zeros(maxconds,2);
itc_results_all=zeros(maxtimepts,maxchans,maxfreqs,maxconds,2);
timepts=zeros(ns,1);

%initialize a bad subject counter(?)
bad_s=zeros(ns,1);

%add prefix if being run on pc
if ispc
    for s_attempt=1:ns
        mgt_coh_list{s_attempt}=[pc_network_drive_prefix,mgt_coh_list{s_attempt}];
    end
end

%load in data
s=0;
for s_attempt=1:ns
    load(mgt_coh_list{s_attempt},'coh_results','mean_data','n_trials','wavelet_evk','wavelet_tot')
    if iscell(coh_results)
        if size(coh_results{1},1) < maxtimepts || size(coh_results{1},3) < maxconds
            bad_s(s_attempt)=1;
            continue
        end
    else
        if size(coh_results,1) < maxtimepts || size(coh_results{1},3) < maxconds
            bad_s(s_attempt)=1;
            continue
        end
    end
    s=s+1;
    if coh_stats_present
        coh_results_all(:,:,:,:,s)=abs(coh_results{1}(1:maxtimepts,:,:,:));
        coh_results_stats(:,:,:,:,s)=coh_results{2}(1:maxtimepts,:,:,:);
    else
        coh_results_all(:,:,:,:,s)=abs(coh_results(1:maxtimepts,:,:,:));
    end
    mean_data_all(:,:,:,s)=mean_data(1:maxtimepts,1:maxchans,:);
    n_trials_all(:,s)=n_trials;
    itc_results_all(:,:,:,:,s)=abs(wavelet_evk(1:maxtimepts,1:maxchans,:,:))./...
        wavelet_tot(1:maxtimepts,1:maxchans,:,:);
end

% baseline-normalize the coherence?


%% prepare dimensional descriptors

%coherence scaling
%coh_abs_min=min(min(min(min(abs(coh_results{1}(:,:,:,:))))));
coh_abs_min=0.5;

%time scaling (in pts for now, can convert from a ms declaration in
%parameters)
time_start_win=1;
time_end_win=maxtimepts;
time_end_max=maxtimepts+2;

%define absolute ms edges
%ms_begin=0;
ms_begin=-500;
ms_end=1500;
%ms_begin=-baseline_size_ms;
%ms_end=max(rearrange_time_ms);
time_ms=linspace(ms_begin,ms_end,time_end_max);

%clip surplus points beyond declared end point
time_ms(maxtimepts+1:end)=[];

%create ms descriptors for the pts
time_xticks=-400:200:1400;
%time_xticks=0:200:1400;
time_xtick=zeros(length(time_xticks),1);
for tick=1:length(time_xticks)
    [c,i]=min(abs(time_ms-time_xticks(tick)));
    time_xtick(tick)=i;
end
[c,time_zero]=min(abs(time_ms));
 
%frequencies %1:4=gamma,5:8=beta,9:12=alpha,13:17=theta,18:20=delta
freqs=param_struct.rate*2./param_struct.scale;
%freqs_range={[18:20],[13:17],[9:12],[5:8],[1:4]};
freqs_range={[19],[14],[10],[6],[3]};
freqs_ytick=[2 6 10 14 18];
freqs_slabel={'Delta','Theta','Alpha','Beta','Gamma'};
freqs_label=round(freqs(freqs_ytick));
freqs_nlabel=cell(maxfreqs,1);
for freq=1:maxfreqs
    freqs_nlabel{freq}=num2str(round(freqs(freq)));
end
freq_colors=distinguishable_colors(length(freqs));

%conditions
cond_inds=param_struct.case_vec;
cond_label_seed={'N50','N10','P10','P50'};
cond_label=cell(maxconds,1);
for cond=1:maxconds
    cond_label{cond}=cond_label_seed{cond_inds(cond)};
end

%subjects
s_label=cell(s,1);
for chosen_s=1:s
    s_label{chosen_s}=num2str(chosen_s);
end
s_colors=distinguishable_colors(s);

%channel pairs
pair_label=cell(maxpairs,1);
for pair=1:maxpairs
    pair_label{pair}=[chan_locs(param_struct.coherence_pairs(pair,1)).labels,'-',chan_locs(param_struct.coherence_pairs(pair,2)).labels];
end
pair_colors=distinguishable_colors(length(pair_label));

%channels
chan_label=cell(maxchans,1);
for chan=1:maxchans
    chan_label{chan}=chan_locs(param_struct.chan_vec(chan)).labels;
end
chan_colors=distinguishable_colors(length(chan_label));

%% data rejection step

s_logic=false(s,3);

%bit to pick out bad data based on high coherence values
% catches channels which are bridged due to high coherence vals, using a
% simple threshold (not well tested)

high_coh_vals=zeros(s,maxpairs);
for chosen_s=1:s
    for pair=1:maxpairs
        chosen_s_data=coh_results_all(:,1:10,:,pair,chosen_s);
        high_coh_vals(chosen_s,pair)=numel(chosen_s_data(chosen_s_data>0.9999)); %/numel(chosen_s_data);
    end
    if sum(high_coh_vals(chosen_s,:),2)>100
        s_logic(chosen_s,1)=false;
    else
        s_logic(chosen_s,1)=true;
    end
end        
high_coh_vals_by_s=sum(high_coh_vals,2);


% bit to pick out bad data based on highly correlated coherence values
% catches channels which are bridged to each other, see subject 3

%figure
%high_coh_corrs=zeros(s,maxpairs);
%coh_rho=zeros(maxpairs,maxpairs,maxfreqs,maxconds);
s_rho=zeros(s,1);
for chosen_s=1:s
    %for freq=14
     %   for cond=1:maxconds
    chosen_s_data=squeeze(coh_results_all(:,10,1,:,chosen_s));
            %chosen_s_data=squeeze(mean(mean(abs(coh_results_all(:,:,:,:,chosen_s)),2),3));
            %coh_rho(:,:,freq,cond)=corr(squeeze(abs(coh_results_all(:,freq,cond,:,chosen_s))));
      %  end
   % end
    rho=corr(chosen_s_data);
    rho=triu(rho,1);
    %s_rho(chosen_s)=numel(rho(rho>.995));
    if numel(rho(rho>.995))>3
        s_logic(chosen_s,2)=false;
    else
        s_logic(chosen_s,2)=true;
    end
end
%imagesc(squeeze(mean(mean(coh_rho,3),4)))
%imagesc(rho)
%caxis([-1 1])
%set(gca,'XTick',1:maxpairs,'XTickLabel',pair_label)
%set(gca,'YTick',1:maxpairs,'YTickLabel',pair_label)
%high_coh_vals(chosen_s,pair)=numel(chosen_s_data(chosen_s_data>0.9999)); %/numel(chosen_s_data);

%bit to exclude subjects which have insufficient data (i.e. not enough
%trials per condition)
trial_thresh=ceil(mean(n_trials_all,2)-std(n_trials_all,0,2));
for chosen_s=1:s
    s_logic(chosen_s,3)=any(n_trials_all(:,chosen_s) > trial_thresh);
end


% combine rejection logic
s_inds=false(s,1);
for chosen_s=1:s
    if all(s_logic(chosen_s,:))
        s_inds(chosen_s)=1;
    end
end

%% bit to plot out bad data

%2,5,17,25, and 32 are examples of high coherence values (>.999)
%3 is an example of highly correlated coherence values

for chosen_s=46 % 3 12 17]
figure(chosen_s)
    for chanpair=1:maxpairs
        plot(coh_results_all(:,9,1,chanpair,chosen_s),'Color',pair_colors(chanpair,:),'LineWidth',2); hold on
    end
hold off
legend(pair_label)
end

%% look at all frequencies for a given condition and pair

chosen_freqs=1:5;
chosen_cond=1;
chosen_pair=14;
chosen_s=1;

figure;
for freq=chosen_freqs
    plot(mean(abs(coh_results_all(:,freqs_range{freq},chosen_cond,chosen_pair,chosen_s)),2),'Color',freq_colors(freq,:)); hold on
    %plot(mean(coh_results_all_test(:,freqs_range{freq},chosen_cond,chosen_pair),2),'Color',freq_colors(freq,:)); hold on
    %plot(mean(coh_results_all(:,freqs_range{freq},chosen_cond,chosen_pair,chosen_s),2),'Color',freq_colors(freq,:)); hold on
end
hold on; plot(ones(maxfreqs,1)*time_zero,linspace(0,1,maxfreqs),'k--'); hold off;
axis([time_start_win time_end_win coh_abs_min 1]);
set(gca,'XTick',time_xtick,'XTickLabel',time_xticks)
legend(freqs_slabel)
title(['Subject ',s_label{chosen_s},' / ',cond_label{chosen_cond},' / ',pair_label{chosen_pair}])

%% look at all freqs for all conditions at one pair

chosen_freqs=1:5;
chosen_pair=14;
chosen_s=1;

figure;
for cond=1:4
    sp(cond)=subplot(2,2,cond);
    for freq=chosen_freqs
        plot(mean(abs(coh_results_all(:,freqs_range{freq},cond,chosen_pair,chosen_s)),2),'Color',freq_colors(freq,:)); hold on
        %plot(mean(coh_results_stats(:,freqs_range{freq},cond,chosen_pair,chosen_s),2),'Color',freq_colors(freq,:)); hold on
        %plot(mean(coh_results_all_test(:,freqs_range{freq},cond,chosen_pair),2),'Color',freq_colors(freq,:)); hold on
    end
    hold off; grid on; axis([time_start_win time_end_win coh_abs_min 1]);
    %hold off; grid on; axis([time_start_win time_end_win 0 16]);
    set(gca,'XTick',time_xtick,'XTickLabel',time_xticks)
    title(['Subject ',s_label{chosen_s},' / ',pair_label{chosen_pair},' / ',cond_label{cond}])
end
linkaxes(sp)
legend(freqs_slabel)

%% look at freqs for all conditions at one pair, for all S's

chosen_freqs=2;
chosen_pair=14;

figure;
for cond=1:maxconds
    sp(cond)=subplot(2,2,cond);
    for freq=chosen_freqs
        for chosen_s=1:s
            %plot(mean(atanh(abs(coh_results_all(:,freqs_range{freq},cond,chosen_pair,chosen_s))),2),'Color',s_colors(chosen_s,:)); hold on
            plot(mean(coh_results_all(:,freqs_range{freq},cond,chosen_pair,chosen_s),2),'Color',s_colors(chosen_s,:)); hold on
            %plot(mean(coh_results_stats(:,freqs_range{freq},cond,chosen_pair,chosen_s),2),'Color',s_colors(chosen_s,:)); hold on
        end
    end
    %hold off; grid on; axis([time_start_win time_end_win .75 4]);
    hold off; grid on; axis([time_start_win time_end_win coh_abs_min 1]);
    set(gca,'XTick',time_xtick,'XTickLabel',time_xticks)
    title([cond_label{cond},' / ',freqs_slabel{chosen_freqs},' / ',pair_label{chosen_pair}])
end
linkaxes(sp)
legend(s_label)

%% look at freqs for all conditions at one pair, averaged across all S's

%delta
%chosen_freqs=18:20;
%theta
chosen_freqs=12:17;
%strongest freq with "that waveform"
%chosen_freqs=16;
%alpha
%chosen_freqs=10:11;
%beta
%chosen_freqs=5:9;
%gamma
%chosen_freqs=1:4;
%all
%chosen_freqs=1:20;

chosen_coh_limits=[0.4 0.8];

for chosen_pair=3
    figure;
    for cond=1:maxconds
        sp(cond)=subplot(2,2,cond);
        for freq=chosen_freqs
            %plot(mean(atanh(abs(coh_results_all(:,freqs_range{freq},cond,chosen_pair,chosen_s))),2),'Color',s_colors(chosen_s,:)); hold on
            %plot(mean(mean(coh_results_all(:,freqs_range{freq},cond,chosen_pair,s_inds),2),5),'Color',freq_colors(freq,:)); hold on; %,'Color',s_colors(chosen_s,:)); hold on
            %calculate data to be plotted (1-d time-series)
            coh_plot_data=nanmean(nanmean(coh_results_all(:,freq,cond,chosen_pair,s_inds),2),5);
            %diff from baseline (4x1 of condition)
            coh_plot_data_base=coh_plot_data - mean(coh_plot_data(1:time_zero));
            %stats
            coh_plot_stats=mean(mean(coh_results_stats(:,freq,cond,chosen_pair,s_inds),2),5)            ;
            %regular
            %plot(coh_plot_data_base,'Color',freq_colors(freq,:)); hold on; %,'Color',s_colors(chosen_s,:)); hold on
            plot(coh_plot_data,'Color',freq_colors(freq,:)); hold on; %,'Color',s_colors(chosen_s,:)); hold on
        end
        %hold off; grid on; axis([time_start_win time_end_win .75 4]);
        grid on; axis([time_start_win time_end_win chosen_coh_limits(1) chosen_coh_limits(2)]);
        plot(ones(10,1)*time_zero,linspace(0,1,10),'k--'); hold off;
        set(gca,'XTick',time_xtick,'XTickLabel',time_xticks)
        title([pair_label{chosen_pair},'/',cond_label{cond},])
    end
linkaxes(sp)
legend(freqs_nlabel(chosen_freqs))
end

%% image coherence in time-freq at a chosen pair

v=zeros(maxpairs,maxconds,2);
for chosen_pair=1:2
    figure(chosen_pair); subplot_dummy=0; 
    for cond=1:maxconds
        subplot_dummy=subplot_dummy+1;
        subplot(2,2,subplot_dummy)
        %imagesc(mean(atanh(coh_results_all(:,:,cond,chosen_pair,s_inds)),5)');
        %imagesc(mean(coh_results_all(:,:,cond,chosen_pair,s_inds),5)');
        contourf(fliplr(nanmean(coh_results_all(:,:,cond,chosen_pair,s_inds),5))');
        shading flat
        %imagesc(mean(coh_results_stats(:,:,cond,chosen_pair,s_inds),5)');
        %axis([time_start_win time_end_win 3 maxfreqs-2]);
        axis([time_start_win time_end_win 1 20]);
        %caxis([1.5 3.5]);
        %caxis([.84 1]);
        %caxis([16 24]);
        v(chosen_pair,subplot_dummy,:) = caxis;
        set(gca,'XTick',time_xtick,'XTickLabel',time_xticks); xlabel('Time (ms)');
        %set(gca,'YTick',freqs_ytick,'YTickLabel',freqs_label); ylabel('Frequency (Hz)');
        set(gca,'YTick',freqs_ytick,'YTickLabel',fliplr(freqs_label)); ylabel('Frequency (Hz)');
        title([pair_label{chosen_pair},' / ',cond_label{cond}])
        hold on; plot(ones(maxfreqs,1)*time_zero,linspace(1,maxfreqs,maxfreqs),'k--'); hold off;
    end
end
v(v==0)=NaN;
cmin=nanmin(nanmin(v(:,:,1))); cmax=nanmax(nanmax(v(:,:,2)));
%cmin=0.45; cmax=1;
for fig=1:2
    figure(fig)
    for splot=1:maxconds; subplot(2,2,splot); caxis([cmin cmax]); end
end

%% plot coherence pairs as lines on a topoplot
data=zeros(maxchans,1);
h=figure;
topoplot(data,chan_locs,'style','blank');
hold on
for pair=1:maxpairs
    line([chan_locs(param_struct.coherence_pairs(pair,1)).topo_x chan_locs(param_struct.coherence_pairs(pair,2)).topo_x],...
        [chan_locs(param_struct.coherence_pairs(pair,1)).topo_y chan_locs(param_struct.coherence_pairs(pair,2)).topo_y],...
        [2.1-randn*.05 2.1+randn*.05],'Color',pair_colors(pair,:)); hold off; %disp(num2str(pair))
end

%% image coherence at a frequency band as a topoplot with colored lines indicating strength

t_start_ms=550;
t_end_ms=650;
f_start_hz=4;
f_end_hz=8;
coh_topo_scale=1;

%convert to points
[c,t_start]=min(abs(time_ms-t_start_ms));
[c,t_end]=min(abs(time_ms-t_end_ms));
[c,f_start]=min(abs(freqs-f_start_hz));
[c,f_end]=min(abs(freqs-f_end_hz));

%define a colormap
cmap=colormap;

figure; subplot_dummy=0; v=zeros(maxconds,2);
dummydata=zeros(length(chan_locs),1);
for cond=1:maxconds
    subplot_dummy=subplot_dummy+1;
    subplot(2,2,subplot_dummy);
    topoplot(dummydata,chan_locs,'style','blank'); hold on;
    for pair=1:maxpairs
        line([chan_locs(param_struct.coherence_pairs(pair,1)).topo_x chan_locs(param_struct.coherence_pairs(pair,2)).topo_x],...
        [chan_locs(param_struct.coherence_pairs(pair,1)).topo_y chan_locs(param_struct.coherence_pairs(pair,2)).topo_y],...
        [2.1-randn*.05 2.1+randn*.05],'Color',...    %coloring happens here
        cmap(round(mean(mean(mean(coh_results_all(t_start:t_end,f_end:f_start,cond,pair,s_inds),1),2),5)*64),:)); hold on;
    end
    hold off; title(sprintf('%s, %d - %d ms, %d - %d Hz',cond_label{cond},t_start_ms(win),t_end_ms(win),f_start_hz,f_end_hz))
end
%cmin=min(v(:,1)); cmax=max(v(:,2));
%cmax=itc_topo_scale;%cmax override
%for splot=1:maxconds; subplot(2,1,splot); caxis([cmin cmax]); end

%% coherence in frequency band as topoplot with colored lines indicating strength, multiple time windows

t_start_ms=[350, 475, 600];
t_end_ms=[350, 475, 600];
f_start_hz=8;
f_end_hz=8;
coh_topo_scale=1;

%convert freqs to points
[c,f_start]=min(abs(freqs-f_start_hz));
[c,f_end]=min(abs(freqs-f_end_hz));

%define a colormap
cmap=colormap;

%define re-scaling constants
linescale=[18,45];

figure('units','normalized','position',[.1 .1 .9 .4]); subplot_dummy=0; %v=zeros(maxconds,2);
dummydata=ones(length(chan_locs),1)*0.3;
coh_linescale_mat=zeros(maxconds,length(t_start_ms),maxpairs);
for cond=1:maxconds
    for win=1:length(t_start_ms)
        [c,t_start]=min(abs(time_ms-t_start_ms(win)));
        [c,t_end]=min(abs(time_ms-t_end_ms(win)));
        subplot_dummy=subplot_dummy+1;
        subplot(maxconds,length(t_start_ms),subplot_dummy)
        topoplot(dummydata,chan_locs,'style','blank','maplimits',[0 1]); hold on;
        for pair=1:maxpairs
            %define [x1 x2], [y1 y2], and [z1 z2] of the arc
            x=[chan_locs(param_struct.coherence_pairs(pair,1)).topo_x chan_locs(param_struct.coherence_pairs(pair,2)).topo_x];
            y=[chan_locs(param_struct.coherence_pairs(pair,1)).topo_y chan_locs(param_struct.coherence_pairs(pair,2)).topo_y];
            %z=[2.1-randn*.05 2.1+randn*.05];
            %determine strength of coherence
            paircoh_for_linecolor=nanmean(nanmean(nanmean(coh_results_all(t_start:t_end,f_end:f_start,cond,pair,s_inds),1),2),5);
            %transform to color
            paircoh_color=paircoh_for_linecolor*64;
            %store color for scaling
            coh_linescale_mat(cond,win,pair)=paircoh_color;
            %scale post hoc
            paircoh_color = round( ( paircoh_color - linescale(1) + 1 ) * 64 / ( linescale(2) - linescale(1) + 1 ) );
            %define its color based on its strength
            linecolor=cmap(paircoh_color,:);
            %direction=mod(pair,2);
            direction=1;
            %plot the arc
            Draw_Arc_Clockwise([x(1),y(1)], [x(2),y(2)], linecolor, 1, direction);
            hold on;
        end
        hold off; title(sprintf('%s, %d - %d ms, %d - %d Hz',cond_label{cond},t_start_ms(win),t_end_ms(win),f_start_hz,f_end_hz))
    end
end
%cmin=min(v(:,1)); cmax=max(v(:,2));
%cmax=itc_topo_scale;%cmax override
%for splot=1:maxconds; subplot(2,1,splot); caxis([cmin cmax]); end
coh_linescale_mat(coh_linescale_mat==0)=NaN;
coh_linescale=[nanmin(nanmin(nanmin(coh_linescale_mat))) nanmax(nanmax(nanmax(coh_linescale_mat)))];

%% image coherence as a topoplot 

%% compare coherence and ITC

%coh_comp_chans=[1 2 4:maxchans];
chosen_chan=9; %8 is CZ, 3 is FZ
chosen_freq=14;
chosen_cond=1;
chosen_pair=4; %14 is FZ-CZ

chan_coh_measure=squeeze(mean(coh_results_all(:,chosen_freq,chosen_cond,chosen_pair,:),1));
chan_coh_stat=squeeze(mean(coh_results_stats(:,chosen_freq,chosen_cond,chosen_pair,:),1));
trial_coh_measure=squeeze(mean(itc_results_all(:,chosen_chan,chosen_freq,chosen_cond,:),1));

figure
scatter3(trial_coh_measure,chan_coh_measure,chan_coh_stat); grid on
xlabel('ITC'); ylabel('Cross-Coherence'); zlabel('Z Score');

%% image ITC in time-freq at a chosen channel

% (16 pair) FZ is 3, FCZ is 16
% (77 pair) FZ is 7, CZ is 16, PZ is 25
chosen_chan=7;

figure; subplot_dummy=0; v=zeros(maxconds,2);
for cond=1:maxconds
    subplot_dummy=subplot_dummy+1;
    subplot(2,2,subplot_dummy)
    itc_plot_data=squeeze(nanmean(itc_results_all(:,chosen_chan,:,cond,s_inds),5))';
    itc_plot_diffmean=itc_plot_data-squeeze(nanmean(nanmean(itc_results_all(:,chosen_chan,:,:,s_inds),4),5))';
    imagesc(itc_plot_diffmean);
    axis([time_start_win time_end_win 1 maxfreqs]);
    v(subplot_dummy,:) = caxis;
    set(gca,'XTick',time_xtick,'XTickLabel',time_xticks); xlabel('Time (ms)');
    set(gca,'YTick',freqs_ytick,'YTickLabel',freqs_label); ylabel('Frequency (Hz)');
    title(['Trial Coherence at ',chan_label{chosen_chan},' : ',cond_label{cond}])
    hold on; plot(ones(maxfreqs,1)*time_zero,linspace(1,maxfreqs,maxfreqs),'k--'); hold off;
end
cmin=min(v(:,1)); cmax=max(v(:,2));
for splot=1:maxconds; subplot(2,2,splot); caxis([cmin cmax]); end

%% image ITC as a scalp plot in a time-frequency window

t_start_ms=0;
t_end_ms=200;
f_start_hz=6;
f_end_hz=10;
itc_topo_scale=1;

%convert to points
[c,t_start]=min(abs(time_ms-t_start_ms));
[c,t_end]=min(abs(time_ms-t_end_ms));
[c,f_start]=min(abs(freqs-f_start_hz));
[c,f_end]=min(abs(freqs-f_end_hz));

figure; subplot_dummy=0; v=zeros(maxconds,2);
for cond=1:maxconds
    subplot_dummy=subplot_dummy+1;
    subplot(2,2,subplot_dummy)
    topo_data=squeeze(mean(mean(mean(itc_results_all(t_start:t_end,:,f_end:f_start,cond,s_inds),1),3),5));
    %topo_data_padded=zeros(length(chan_locs),1);
    %for chan=1:length(good_inds)
    %    topo_data_padded(good_inds(chan))=topo_data(chan);
    %end
    topoplot(topo_data,chan_locs,'maplimits',[0 itc_topo_scale]);
    %v(subplot_dummy,:) = caxis;
    title(sprintf('%d - %d ms, %d - %d Hz',t_start_ms,t_end_ms,f_start_hz,f_end_hz))
end
%cmin=min(v(:,1)); cmax=max(v(:,2));
%cmax=itc_topo_scale;%cmax override
%for splot=1:maxconds; subplot(2,1,splot); caxis([cmin cmax]); end

%% image ITC as a scalp plot in a series of time-frequency windows

t_start_ms=[-200, 0, 200, 400, 600];
t_end_ms=[0, 200, 400, 600, 800];
f_start_hz=4;
f_end_hz=8;
itc_topo_scale=0.6;

figure; subplot_dummy=0; v=zeros(maxconds,2);
%convert freqs to pts
[c,f_start]=min(abs(freqs-f_start_hz));
[c,f_end]=min(abs(freqs-f_end_hz));
for cond=1:maxconds
    for win=1:length(t_start_ms)
        %convert times to points
        [c,t_start]=min(abs(time_ms-t_start_ms(win)));
        [c,t_end]=min(abs(time_ms-t_end_ms(win)));
        %plot
        subplot_dummy=subplot_dummy+1;
        subplot(maxconds,length(t_start_ms),subplot_dummy)
        topo_data=squeeze(nanmean(nanmean(nanmean(itc_results_all(t_start:t_end,:,f_end:f_start,cond,s_inds),1),3),5));
        %topo_data_padded=zeros(length(chan_locs),1);
        %for chan=1:length(good_inds)
        %    topo_data_padded(good_inds(chan))=topo_data(chan);
        %end
        topoplot(topo_data,chan_locs,'maplimits',[0 itc_topo_scale]);
        %v(subplot_dummy,:) = caxis;
        title(sprintf('%s, %d - %d ms, %d - %d Hz',cond_label{cond},t_start_ms(win),t_end_ms(win),f_start_hz,f_end_hz))
    end
end
%cmin=min(v(:,1)); cmax=max(v(:,2));
%cmax=itc_topo_scale;%cmax override
%for splot=1:subplot_dummy; subplot(maxconds,length(t_start_ms),splot); caxis([cmin cmax]); end