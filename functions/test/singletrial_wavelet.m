%% set params

chan = 7; %FZ
pair = 57;
pair_chans = opt.coherence_pairs(pair,:);
freq = 12; %theta

cond=1;
dotted_times = [34,113];

srate = 256; lp_cutoff = 16; %filter settings

n_conds = size(output.trials,2);
conds = {'Gain', 'Loss'};

cmap = repmat(linspace(0.25,1,256)',[1 3]); %"grayscale map" for base map
cmap = circularize_map(cmap);

%% verbose time scaling section

ds_ratio = 2;
n_timepts = size(dataWG{1},1);

ms_start=-opt.prestim_ms;
ms_end=ms_start + ( n_timepts * (1000 / srate*ds_ratio) );

t_ms = linspace(ms_start,ms_end,n_timepts+1);
t_ms = t_ms(1:end-1);

t_start_b_ms=-500; t_end_b_ms=-200;
[~,t_start_b]=min(abs(t_ms-t_start_b_ms));
[~,t_end_b]=min(abs(t_ms-t_end_b_ms));

ms_tickint=200;
t_xtick_ms=(ms_start+mod(ms_start,ms_tickint)): ...
    ms_tickint:(ms_end-mod(ms_end,ms_tickint));
t_xtick=zeros(length(t_xtick_ms),1);
for tick=1:length(t_xtick_ms)
    [~,i]=min(abs(t_ms-t_xtick_ms(tick)));
    t_xtick(tick)=i;
end
[~,t_zero]=min(abs(t_ms));
t_xtick_ms = t_xtick_ms ./ 1000;

%% put the wavelet TFR into matrices

amp = cell(1, n_conds); ampnorm = cell(1, n_conds);
amp1 = cell(1, n_conds); amp2 = cell(1, n_conds);
pow = cell(1, n_conds); phi = cell(1, n_conds);
phi_comp = cell(1, n_conds); 
phi1 = cell(1, n_conds); phi2 = cell(1, n_conds);
phi_all = [];
for cond = 1:n_conds
    data = squeeze(dataWG{cond}(:,:,chan,freq));
    amp{cond} = real(data);
    ampnorm{cond} = bsxfun(@rdivide, real(data), max(abs(real(data))) );
    pow{cond} = data.*conj(data);
    phi{cond} = angle(data);
    phi_comp{cond} = data;
    amp1{cond} = real(squeeze(dataWG{cond}(:,:,pair_chans(1),freq)));
    amp2{cond} = real(squeeze(dataWG{cond}(:,:,pair_chans(2),freq)));
    phi1{cond} = angle(squeeze(dataWG{cond}(:,:,pair_chans(1),freq)));
    phi2{cond} = angle(squeeze(dataWG{cond}(:,:,pair_chans(2),freq)));
    phi_all = [phi_all, angle(data)];
end


%% filter erptrial data and assign to a new matrix

erptrial_filt = permute(output.erptrial,[1 3 2]);
erptrial_filt = filter_erpdata(erptrial_filt, srate, lp_cutoff);
erptrial_filt = permute(erptrial_filt,[1 3 2]);

dataE = cell(1, n_conds);
for cond = 1:n_conds
    dataE{cond} = erptrial_filt(:,output.trials(:,cond),chan);
end

%% quickly check up on the ERO vals

figure;
for cond=1:2
    subplot(1,2,cond);
    itcplotdata = flipud(10*log10( bsxfun(@rdivide, squeeze(output.wave_tot(:,chan,:,cond)).^2, ...
        mean(squeeze(output.wave_tot(t_start_b:t_end_b,chan,:,cond)).^2,1)) )');
    contourf( itcplotdata );
end

%% quickly check up on the ITC vals

figure;
for cond=1:2
    subplot(1,2,cond);
    itcplotdata = flipud(abs(squeeze(output.wave_evknorm(:,chan,:,cond)))');
    contourf( itcplotdata );
end

%% quickly check up on the ISPC vals

figure;
for cond=1:2
    subplot(1,2,cond);
    ispcplotdata = flipud(abs(squeeze(output.coh(:,:,cond,pair)))');
    contourf( ispcplotdata );
end

%% plot the trial data as a butterfly plot

figure;
subplot(211);
plot(dataE{1});
subplot(212);
plot(dataE{2});

%% plot the bandpass filtered trials

figure;
subplot(211);
plot(amp{1});
subplot(212);
plot(amp{2});

%% plot the bandpass filtered trials
% (inter-site, separately)

figure;
subplot(211);
plot(amp1{1}+2,'k'); hold on;
plot(amp2{1}-2,'m');
subplot(212);
plot(amp1{2}+2,'k'); hold on;
plot(amp2{2}-2,'m');

%% plot the bandpass filtered trials
% (inter-site, as difference)

figure;
subplot(211);
plot(amp1{1} - amp2{1});
subplot(212);
plot(amp1{2} - amp2{2});

%% plot "normalized" bandpass filtered trials

figure;
subplot(211);
plot(ampnorm{1}); ylim([-2 2]);
subplot(212);
plot(ampnorm{2}); ylim([-2 2]);

%% plot a "stacked" visualization of the bandpass-filtered trials

figure
for cond=1:2
n_trials = size(amp{cond},2);

subplot(1,2,cond);
for trial=1:n_trials
    plot(amp{cond}(:,trial)+trial);
    hold on;
end
plot(ones(n_timepts,1)*dotted_times(1),linspace(0,n_trials,n_timepts),'k--');
plot(ones(n_timepts,1)*dotted_times(2),linspace(0,n_trials,n_timepts),'k--');
hold off
ylim([-1 n_trials+2]);
end

%% plot a "stacked" visualization of the bandpass-filtered trials
% (inter-site)

figure;
hold on;
trial_colors = distinguishable_colors(n_trials);
for trial=1:n_trials
    plot(amp1{cond}(:,trial)+trial,'Color',trial_colors(trial,:));
    plot(amp2{cond}(:,trial)+trial,'Color',trial_colors(trial,:));
end
plot(ones(n_timepts,1)*dotted_times(1),linspace(0,n_trials,n_timepts),'k--');
plot(ones(n_timepts,1)*dotted_times(2),linspace(0,n_trials,n_timepts),'k--');
hold off
ylim([-1 n_trials+2]);

%% plot a "phase image" where trials is the y axis and time is the x axis

figure;
for cond=1:2
    subplot(1,2,cond)
    imagesc(phi{cond}')
    %contourf(phi{cond}');
    colormap(cmap);
    c=colorbar('YLim',[-pi pi],'YTick',[-pi,0,pi],'YTickLabel',{'-\pi','0','\pi'}, ...
        'FontSize', 12);
    hold on;
    plot(ones(n_timepts,1)*dotted_times(1),linspace(0,n_trials,n_timepts),'k--');
    plot(ones(n_timepts,1)*dotted_times(2),linspace(0,n_trials,n_timepts),'k--');
    hold off;
    set(gca,'XTick',t_xtick,'XTickLabel',t_xtick_ms); xlabel('Time (ms)');
    ylabel('Trial #');
    title(conds{cond});
end
axis_prunelabels;

%% (inter-site)

figure;
for cond=1:2
    subplot(1,2,cond)
    imagesc( deunwrap( phi1{cond} - phi2{cond} )' )
    %contourf(phi{cond}');
    colormap(cmap);
    c=colorbar('YLim',[-pi pi],'YTick',[-pi,0,pi],'YTickLabel',{'-\pi','0','\pi'}, ...
        'FontSize', 12);
    hold on;
    plot(ones(n_timepts,1)*dotted_times(1),linspace(0,n_trials,n_timepts),'k--');
    plot(ones(n_timepts,1)*dotted_times(2),linspace(0,n_trials,n_timepts),'k--');
    hold off;
    set(gca,'XTick',t_xtick,'XTickLabel',t_xtick_ms); xlabel('Time (ms)');
    ylabel('Trial #');
    title(conds{cond});
end
axis_prunelabels;

%% interleaved phase image

trials_terse = output.trials;
trials_terse( sum(trials_terse,2) == 0, :) = [];
map = [2 1];

figure;
% phasemaps
for cond=1:2
    subplot(1,2,cond)
    phasemap = phi_all;
    phasemap(:, trials_terse(:,map(cond)) ) = 0;
    imagesc(phasemap');
    set(gca, 'YDir', 'normal');
    %contourf(phi{cond}');
    colormap(cmap);
    c=colorbar('YLim',[-pi pi],'YTick',[-pi,0,pi],'YTickLabel',{'-\pi','0','\pi'}, ...
        'FontSize', 12);
    hold on;
    plot(ones(n_timepts,1)*dotted_times(1),linspace(1,length(trials_terse),n_timepts),'y--');
    plot(ones(n_timepts,1)*dotted_times(2),linspace(1,length(trials_terse),n_timepts),'m--');
    plot(ones(n_timepts,1)*t_zero,linspace(1,length(trials_terse),n_timepts),'w');
    hold off;
    set(gca,'XTick',t_xtick,'XTickLabel',t_xtick_ms); xlabel('Time (s)');
    ylabel('Trial #');
    title(conds{cond});
end

%% ITC over time to correspond with above

figure;
for cond=1:2
    subplot(1,2,cond);
    plot(abs(output.wave_evknorm(:,chan,freq,cond)))
    set(gca,'XTick',t_xtick,'XTickLabel',t_xtick_ms); xlabel('Time (s)');
    ylabel('ITC')
    ylim([0 0.8])
end
set_print_size(16,2);

%% interleaved phase image
% (inter-site)

trials_terse = output.trials;
trials_terse( sum(trials_terse,2) == 0, :) = [];
map = [2 1];

figure;
for cond=1:2
    subplot(1,2,cond)
    phasemap = phi_all;
    phasemap(:, trials_terse(:,map(cond)) ) = 0;
    imagesc(phasemap');
    set(gca, 'YDir', 'normal');
    %contourf(phi{cond}');
    colormap(cmap);
    c=colorbar('YLim',[-pi pi],'YTick',[-pi,0,pi],'YTickLabel',{'-\pi','0','\pi'}, ...
        'FontSize', 12);
    hold on;
    plot(ones(n_timepts,1)*dotted_times(1),linspace(1,length(trials_terse),n_timepts),'y--');
    plot(ones(n_timepts,1)*dotted_times(2),linspace(1,length(trials_terse),n_timepts),'m--');
    plot(ones(n_timepts,1)*t_zero,linspace(1,length(trials_terse),n_timepts),'w');
    hold off;
    set(gca,'XTick',t_xtick,'XTickLabel',t_xtick_ms); xlabel('Time (s)');
    ylabel('Trial #');
    title(conds{cond});
end

%% plot a polar plot of trial phases at a given time point

for cond = 1:2
    figure;
    for time = 1:length(dotted_times)
        subplot(1,2,time);
        polar([phi{cond}(dotted_times(time),:); phi{cond}(dotted_times(time),:)], ...
            [zeros(size(phi{cond}(dotted_times(time),:))); ones(size(phi{cond}(dotted_times(time),:)))], ...
            'k-');
        itcval = abs(output.wave_evknorm(dotted_times(time),chan,freq,cond));
        title(['ITC = ',num2str( itcval )]);
    end
end
%% (inter site)

for cond = 1:2
    figure;
    for time = 1:length(dotted_times)
        subplot(1,2,time);
        phase_diffs = phi1{cond}(dotted_times(time),:) - phi2{cond}(dotted_times(time),:);
        polar([phase_diffs; phase_diffs], ...
            [zeros(size(phase_diffs)); ones(size(phase_diffs))], ...
            'k-');
        ispcval = abs(output.coh(dotted_times(time),freq,cond,pair));
        title(['ISPC = ',num2str( ispcval )]);
    end
end


%% plot a histogram of trial phases at a given time point

sp_rowlabel = conds;
sp_columnlabel = { [num2str(t_ms(dotted_times(1)),3),' ms'], [num2str(t_ms(dotted_times(2)),3),' ms'] };
sp_xlabel = 'Time Region';
sp_ylabel = 'Condition';
overtitle = '';
sp_dims = [2 2];

raxis_lim = 0.6;

pts = linspace(-3.2, 3.2, 100);

spd=0;
figure;
for cond = 1:2
    for time = 1:length(dotted_times)
        
        spd=spd+1;
        subplot(2,2,spd);
        
        %cheat to make the raxis big enough
        t = 0 : .01 : 2 * pi;
        P = polar2(t, raxis_lim * ones(size(t)));
        set(P, 'Visible', 'off')
        hold on;
        
        %histogram( phi{cond}(dotted_times(time),:), 10 );
        %axis( [-4 4 0 25] );
        
        [f,xi] = ksdensity(phi{cond}(dotted_times(time),:), pts);
        %plot(xi,f);
        %axis( [-pi pi 0 0.6] );
        
        polar( xi, f );
        
        %itcval = abs(output.wave_evknorm(dotted_times(time),chan,freq,cond));
        %text(-1.5*raxis_lim,.75*raxis_lim,['ITC = ',num2str( itcval, 3 )]);
        
    end
end
%adorn_plots(sp_rowlabel, sp_columnlabel, sp_xlabel, sp_ylabel, overtitle, sp_dims);

%% plot a compass plot of trial phases at a given time point

for cond = 1:2
    figure;
    for time = 1:length(dotted_times)
        subplot(1,2,time);
        compass( phi_comp{cond}(dotted_times(time),:) );
        itcval = abs(output.wave_evknorm(dotted_times(time),chan,freq,cond));
        title(['ITC = ',num2str( itcval )]);
    end
end

%% plot a rose of trial phases at a given time point

sp_rowlabel = conds;
sp_columnlabel = { [num2str(t_ms(dotted_times(1)),3),' ms'], [num2str(t_ms(dotted_times(2)),3),' ms'] };
sp_xlabel = 'Time Region';
sp_ylabel = 'Condition';
overtitle = '';
sp_dims = [2 2];

raxis_lim = 12;

spd=0;
figure;
for cond = 1:2
    for time = 1:length(dotted_times)
        
        spd=spd+1;
        subplot(2,2,spd);
        
        %cheat to make the raxis big enough
        t = 0 : .01 : 2 * pi;
        P = polar2(t, raxis_lim * ones(size(t)));
        set(P, 'Visible', 'off')
        hold on;
        
        %rose
        rose2( phi{cond}(dotted_times(time),:) );
        %hold off;
        
        itcval = abs(output.wave_evknorm(dotted_times(time),chan,freq,cond));
        text(-1.5*raxis_lim,.75*raxis_lim,['ITC = ',num2str( itcval, 3 )]);
        
    end
end
adorn_plots(sp_rowlabel, sp_columnlabel, sp_xlabel, sp_ylabel, overtitle, sp_dims);

%% plot a rose of trial phases at a given time point
% (inter-site)

raxis_lim = 15;

for cond = 1:2
    figure;
    
    for time = 1:length(dotted_times)
        subplot(1,2,time);
        
        phase_diffs = phi1{cond}(dotted_times(time),:) - phi2{cond}(dotted_times(time),:);
        
        %cheat to make the raxis big enough
        t = 0 : .01 : 2 * pi;
        P = polar(t, raxis_lim * ones(size(t)));
        set(P, 'Visible', 'off')
        hold on;
        
        %rose
        rose2( phase_diffs );
        %hold off;
        
        ispcval = abs(output.coh(dotted_times(time),freq,cond,pair));
        title(['ISPC = ',num2str( ispcval )]);
        
    end
    plottitle(conds{cond},2)
end

%% plot a wind rose where the levels are conditions??

di_bins = [.9, 1.5, 2.1];

windshape = [ phi{1} phi{2} ]*180/pi;
wind_di = [ ones(size(phi{1})) ones(size(phi{2}))*2 ];

% red = Loss, blue = gain
for time = 1:length(dotted_times)
    figure;
    
    wind_rose(windshape(dotted_times(time),:), wind_di(dotted_times(time),:), ...
        'n', 18, 'ci', [10 20], 'lablegend', 'condition', 'quad', 3, ...
        'di', di_bins);
    
    itcval = abs(output.wave_evknorm(dotted_times(time),chan,freq,cond));
    title([num2str(t_ms(dotted_times(time)),3),'ms , ITC = ',num2str( itcval )]);

end

%% plot a wind rose where the levels are times???

di_bins = [.9, 1.5, 2.1];

% red = 350 ms, blue = -250 ms
for cond = 1:n_conds
    figure;
    
    windshape = [ phi{cond}(dotted_times(1),:) phi{cond}(dotted_times(2),:) ]*180/pi;
    wind_di = [ ones(size(phi{cond}(dotted_times(1),:))) ones(size(phi{cond}(dotted_times(2),:)))*2 ];
    
    wind_rose(windshape, wind_di, ...
        'n', 18, 'ci', [10 20], 'lablegend', 'condition', 'quad', 3, ...
        'di', di_bins);
    
    itcval = abs(output.wave_evknorm(dotted_times(time),chan,freq,cond));
    title([conds{cond},', ITC = ',num2str( itcval )]);

end

%% plot the power, as single-trial butterflies and as trial averages

powlog = 10 * log10( bsxfun(@rdivide, pow, meanx(pow(t_start_b:t_end_b, :), 2)) );
powlogmean = 10 * log10( bsxfun(@rdivide, mean(pow,2), meanx(pow(t_start_b:t_end_b, :), [])) );

figure;
subplot(2,2,1);
plot(pow);
subplot(2,2,2);
plot(powlog);
subplot(2,2,3);
plot(mean(pow,2));
subplot(2,2,4);
plot(powlogmean);

%%

figure;
plot(itc);