%% pre-load procedure (need only be run once)

%load('/export/home/mike/matlab/batch/ac252_eec_325/ac252_eec_ac252_eec_opt.mat')
load('/export/home/mike/matlab/batch/ac252_eec_fpo90/ac252_eec_fpo90_opt.mat')
load('/export/home/mike/matlab/origin/coords/61chans_ns.mat')

%mat_dir = '/active_projects/mike/ac252_eec/';
mat_dir = '/active_projects/mike/ac252_eec_fpo90/';

%mat_list = strsplit(ls(mat_dir),'\n'); mat_list = mat_list(1:end-1);
dir_contents = dir(mat_dir);
mat_list_r = {dir_contents(3:end).name}';
n_files = length(mat_list_r);
uniq_id = cell(n_files, 1);
for f = 1:n_files
    uniq_id{f} = mat_list_r{f}(7:17);
end
s_demogs_r = table(uniq_id, mat_list_r);
s_demogs2 = s_demogs(s_inds_g(:,1),:);
[a,b] = join(s_demogs2, s_demogs_r, 'keys', 1);
mat_list_r = a.mat_list_r;
s_inds_g_r = s_inds_g(s_inds_g(:,1),:);


%% load

freqs = 2*opt.rate ./ opt.wavelet_scales;
freqs = round(freqs*10)/10;
freqsf = fliplr(freqs);
n_freqs = length(freqs);

n_chans = length(opt.chan_vec);
n_pairs = size(opt.coherence_pairs,1);

rcohdata = zeros(n_pairs, n_freqs, n_files);
rpowdata = zeros(n_chans, n_freqs, n_files);

for f=1:n_files
    
    load(fullfile(mat_dir,mat_list_r{f}));
    
    rcohdata(:,:,f) = abs(coh);
    rpowdata(:,:,f) = wave_totpow;
    
end

clear mat_dir dir_contents uniq_id freqs freqsf n_freqs n_chans ...
    n_pairs f n_files n_interpchans n_trials trials coh wave_totpow a b
    

%%

figure
for f=1:n_files
    
    plot( fliplr( squeeze(rpowdata(:,:,f)) )');
    set(gca,'xtick',1:2:imp.maxfreqs,'xticklabel',freqsf(1:2:imp.maxfreqs));
    pause;
    clf;
    
end

%%

figure
for f=1:n_files
    
    plot( fliplr( squeeze(rcohdata(:,:,f)) )');
    axis([1 31 0 1]);
    pause;
    clf;
    
end

%% power spectra

figure
%plot( fliplr( squeeze(mean(powdata, 3)) )');
%powspectra_data = 10.*log10( bsxfun(@rdivide, powdata, freqs));
%powspectra_data = 10.*log10( powdata );
powspectra_data = rpowdata;
%powspectra_data = log(rpowdata);

plot( fliplr( squeeze(mean(powspectra_data, 3)) )');
set(gca,'xtick',1:imp.maxfreqs,'xticklabel',flipud(scl.f_nlabel));
xlabel('Frequency (Hz)');
ylabel('Power')
%ylabel('log(Power)')


%% topography of power spectra

%topolimits = [0 4];
%topolimits = [-16 3];

%cmap = parula(256);
%cmap = makecmap(topolimits, Inf);
cmap = pmkmp(256, 'cubicl');

sp_d = numSubplots(imp.maxfreqs);
figure;
v = zeros(imp.maxfreqs, 2);
for freq_range=1:4
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    %[~,f_ind]=min(abs(scl.freqs-pp.f_indiv_hz(f)));
    subplot(2, 2, freq_range);
    %powtopo_data = squeeze(meanx( 10.*log10(rpowdata(:,f_end:f_start,:)), 1));
    %powtopo_data = squeeze(meanx(rpowdata(:,f_end:f_start,:), 1));
    powtopo_data = squeeze(meanx(log(rpowdata(:,f_end:f_start,:)), 1));
    v(freq_range, 1) = min(powtopo_data); v(freq_range, 2) = max(powtopo_data);
    h = topoplot( powtopo_data, chan_locs, ...
        'colormap', cmap, 'style', 'fill', 'numcontour', 8, ...
        'electrodes', 'off', 'maplimits', 'maxmin');
    set(h,'EdgeColor','None');
    colorbar;
    subtitle = sprintf('%1.1f - %1.1f Hz',scl.freqs(f_start), ...
        scl.freqs(f_end) );
    title(subtitle);
end
%c(1) = min(v(:,1)); c(2) = max(v(:,2));
%figure;
%colorscale_plot(topolimits, cmap, 0.5);

%% intra vs. inter-regional pairs

p_set={1:41, 42:90};

figure
for s=1:length(p_set)

subplot(1,2,s);
rcoh_plotdata = fliplr( squeeze(mean(rcohdata(p_set{s},:,:), 3)) )';
plot(rcoh_plotdata);
axis([1 imp.maxfreqs 0 0.7]); grid on;
set(gca, 'xtick', 1:imp.maxfreqs, 'xticklabel', flipud(scl.f_nlabel));
xlabel('Frequency (Hz'); ylabel('Coherence');
clickableLegend(scl.p_label(p_set{s}));
end

%% regional hypotheses

figure
for hyp=1:length(opt.pair_indlbls)
subplot(2,3,hyp);
pair_hypinds = opt.pair_inds == hyp;
rcoh_plotdata = fliplr( squeeze(mean(rcohdata(pair_hypinds,:,:), 3)) )';
plot(rcoh_plotdata);
axis([1 imp.maxfreqs 0 0.7]); grid on;
set(gca, 'xtick', 1:imp.maxfreqs, 'xticklabel', flipud(scl.f_nlabel));
xlabel('Frequency (Hz'); ylabel('Coherence');
clickableLegend(scl.p_label(pair_hypinds));
title(opt.pair_indlbls{hyp});
end

%% scatter resting intra-regional power vs. coherence

axis_lims = [0 1 -4 4];
%axis_lims = [0 1 -4 4];

for hyp=4:6
figure;
pair_hypinds = opt.pair_inds == hyp;
pair_hypchans = unique(opt.coherence_pairs(pair_hypinds,:));
overtitle=opt.pair_indlbls{hyp};
for freq_range=1:4
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    [~,f_ind]=min(abs(scl.freqs-pp.f_indiv_hz(freq_range)));
    subplot(1,4,freq_range);
    subtitle = sprintf('%1.1f - %1.1f Hz',scl.freqs(f_start), ...
        scl.freqs(f_end) );
    hold on;
for group=pp.chosen_g(pp.plotn_g)
    
    x_scatterdata = squeeze(meanx(rcohdata(pair_hypinds, f_end:f_start, s_inds_g_r(:,group)),3));
    %y_scatterdata = squeeze(meanx(rpowdata(pair_hypchans, f_end:f_start, s_inds_g_r(:,group)),3));
    y_scatterdata = squeeze(meanx(log(rpowdata(pair_hypchans, f_end:f_start, s_inds_g_r(:,group))),3));
    
    plot( x_scatterdata, y_scatterdata, [scl.g_color{group},'.'] );
    plot( mean(x_scatterdata), mean(y_scatterdata), [scl.g_color{group},'x'], ...
        'markersize', 15, 'linewidth', 4);
    
end
axis(axis_lims); grid on;
xlabel('Coherence'); ylabel('log(Power)');
title(subtitle);
plot(linspace(axis_lims(1),axis_lims(2),100),linspace(axis_lims(3),axis_lims(4),100),'k--');
hold off;
end
%tightfig;
plottitle(overtitle,2);
set(gcf, 'Position', [1300 100 1900 325]);
end

%% scatter resting and baseline coherence

baseline_region=[-500 -200];
[~,t_start_b]=min(abs(scl.t_ms-baseline_region(1)));
[~,t_end_b]=min(abs(scl.t_ms-baseline_region(2)));

for hyp=1
figure;
pair_hypinds = opt.pair_inds == hyp;
overtitle=opt.pair_indlbls{hyp};
for freq_range=1:4
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    [~,f_ind]=min(abs(scl.freqs-pp.f_indiv_hz(freq_range)));
    subplot(1,4,freq_range);
    subtitle = sprintf('%1.1f - %1.1f Hz',scl.freqs(f_start), ...
        scl.freqs(f_end) );
    hold on;
for group=pp.chosen_g(pp.plotn_g)
    
    x_scatterdata = squeeze(meanx(rcohdata(pair_hypinds, f_end:f_start, s_inds_g_r(:,group)),3));
    y_scatterdata = meanx( cohdata(t_start_b:t_end_b, f_end:f_start, :, pair_hypinds, s_inds_g(:,group)), 5);
    %x_scatterdata = squeeze(meanx(rcohdata(pair_hypinds, f_ind, s_inds_g_r(:,group)),3));
    %y_scatterdata = meanx( cohdata(t_start_b:t_end_b, f_ind, :, pair_hypinds, s_inds_g(:,group)), 5);
    
    scatter( x_scatterdata, y_scatterdata, scl.g_color{group} );
    
end
axis([0 1 0 1]);
xlabel('Resting'); ylabel('Baseline');
title(subtitle);
plot(linspace(0,1,100),linspace(0,1,100),'k--');
hold off;
end
%tightfig;
plottitle(overtitle,2);
set(gcf, 'Position', [1300 100 1900 325]);
end

%% scatter baseline and event-related coherence

baseline_region=[-500 -200];
[~,t_start_b]=min(abs(scl.t_ms-baseline_region(1)));
[~,t_end_b]=min(abs(scl.t_ms-baseline_region(2)));

event_region=[200 400];
[~,t_start]=min(abs(scl.t_ms-event_region(1)));
[~,t_end]=min(abs(scl.t_ms-event_region(2)));

for hyp=1
figure;
pair_hypinds = opt.pair_inds == hyp;
overtitle=opt.pair_indlbls{hyp};
for freq_range=1:4
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    [~,f_ind]=min(abs(scl.freqs-pp.f_indiv_hz(freq_range)));
    subplot(1,4,freq_range);
    subtitle = sprintf('%1.1f - %1.1f Hz',scl.freqs(f_start), ...
        scl.freqs(f_end) );
    hold on;
for group=pp.chosen_g(pp.plotn_g)
    
    x_scatterdata = meanx( cohdata(t_start_b:t_end_b, f_end:f_start, :, pair_hypinds, s_inds_g(:,group)), 5);
    y_scatterdata = meanx( cohdata(t_start:t_end, f_end:f_start, :, pair_hypinds, s_inds_g(:,group)), 5);
    %x_scatterdata = meanx( cohdata(t_start_b:t_end_b, f_ind, :, pair_hypinds, s_inds_g(:,group)), 5);
    %y_scatterdata = meanx( cohdata(t_start:t_end, f_ind, :, pair_hypinds, s_inds_g(:,group)), 5);
    
    scatter( x_scatterdata, y_scatterdata, scl.g_color{group} );
    
end
axis([0 1 0 1]);
xlabel('Baseline'); ylabel('Event');
title(subtitle);
plot(linspace(0,1,100),linspace(0,1,100),'k--');
hold off;
end
%tightfig;
plottitle(overtitle,2);
set(gcf, 'Position', [1300 100 1900 325]);
end

%% scatter baseline and event-related coherence

baseline_region=[-500 -200];
[~,t_start_b]=min(abs(scl.t_ms-baseline_region(1)));
[~,t_end_b]=min(abs(scl.t_ms-baseline_region(2)));

event_region=[200 400];
[~,t_start]=min(abs(scl.t_ms-event_region(1)));
[~,t_end]=min(abs(scl.t_ms-event_region(2)));

for hyp=1
figure;
pair_hypinds = opt.pair_inds == hyp;
overtitle=opt.pair_indlbls{hyp};
for freq_range=1:4
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    [~,f_ind]=min(abs(scl.freqs-pp.f_indiv_hz(freq_range)));
    subplot(1,4,freq_range);
    subtitle = sprintf('%1.1f - %1.1f Hz',scl.freqs(f_start), ...
        scl.freqs(f_end) );
    hold on;
for group=pp.chosen_g(pp.plotn_g)
    
    x_scatterdata = meanx( cohdata(t_start_b:t_end_b, f_end:f_start, :, pair_hypinds, s_inds_g(:,group)), 5);
    y_scatterdata = meanx( cohdata(t_start:t_end, f_end:f_start, :, pair_hypinds, s_inds_g(:,group)), 5);
    %x_scatterdata = meanx( cohdata(t_start_b:t_end_b, f_ind, :, pair_hypinds, s_inds_g(:,group)), 5);
    %y_scatterdata = meanx( cohdata(t_start:t_end, f_ind, :, pair_hypinds, s_inds_g(:,group)), 5);
    
    scatter( x_scatterdata, y_scatterdata, scl.g_color{group} );
    
end
axis([0 1 0 1]);
xlabel('Baseline'); ylabel('Event');
title(subtitle);
plot(linspace(0,1,100),linspace(0,1,100),'k--');
hold off;
end
%tightfig;
plottitle(overtitle,2);
set(gcf, 'Position', [1300 100 1900 325]);
end

%% scatter resting, baseline, and event-related coherence

axis_lims=[0 1];
axis_lims_z=[0 1];
%axis_lims_z=[-.3 .3];

for hyp=1:6
figure;
pair_hypinds = opt.pair_inds == hyp;
overtitle=opt.pair_indlbls{hyp};
for freq_range=1:4
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    %[~,f_ind]=min(abs(scl.freqs-pp.f_indiv_hz(freq_range)));
    subplot(1,4,freq_range);
    subtitle = sprintf('%1.1f - %1.1f Hz',scl.freqs(f_start), ...
        scl.freqs(f_end) );
    hold on;
for group=pp.chosen_g(pp.plotn_g)
    
    % rest
    x_scatterdata = squeeze(meanx(rcohdata(pair_hypinds, f_end:f_start, s_inds_g_r(:,group)),3));
    % baseline
    y_scatterdata = meanx( cohdata(t_start_b:t_end_b, f_end:f_start, :, pair_hypinds, s_inds_g(:,group)), 5);
    % event
    z_scatterdata = meanx( cohdata(t_start:t_end, f_end:f_start, :, pair_hypinds, s_inds_g(:,group)), 5);
    %z_scatterdata = meanx( cohdata(t_start:t_end, f_end:f_start, pp.cond_diff{1}, pair_hypinds, s_inds_g(:,group)), 5) - ...
    %    meanx( cohdata(t_start:t_end, f_end:f_start, pp.cond_diff{2}, pair_hypinds, s_inds_g(:,group)), 5);
    %z_scatterdata = meanx( cohdata(t_start:t_end, f_end:f_start, 2, pair_hypinds, s_inds_g(:,group)), 5) - ...
    %    y_scatterdata;
    
    %x_scatterdata = squeeze(meanx(rcohdata(pair_hypinds, f_ind, s_inds_g_r(:,group)),3));
    %y_scatterdata = meanx( cohdata(t_start_b:t_end_b, f_ind, :, pair_hypinds, s_inds_g(:,group)), 5);
    
    scatter3( x_scatterdata, y_scatterdata, z_scatterdata, [scl.g_color{group},'.'] );
    plot3( mean(x_scatterdata), mean(y_scatterdata), mean(z_scatterdata), ...
        [scl.g_color{group},'x'], 'markersize', 15, 'linewidth', 4);
    
    [x,y,z] = fitline_3dpts(x_scatterdata, y_scatterdata, z_scatterdata, [-2 2]);
    plot3(x,y,z,scl.g_color{group});
    
end
axis([axis_lims(1) axis_lims(2) axis_lims(1) axis_lims(2) axis_lims_z(1) axis_lims_z(2)]); grid on;
xlabel('Resting'); ylabel('Baseline'); zlabel('Event');
title(subtitle);
plot3(linspace(axis_lims(1),axis_lims(2),100),linspace(axis_lims(1),axis_lims(2),100), ...
    linspace(axis_lims_z(1),axis_lims_z(2),100),'k--');
hold off;
%view(-20, 25);
view(0, 90); %resting = x, baseline = y
%view(90, 0); %baseline = x, event = y
%view(0, 0); % resting = x, event = y
end
%tightfig;
plottitle(overtitle,2);
set(gcf, 'Position', [1300 100 1900 325]);
end

%% scatter resting, baseline, and event-related POWER

baseline_region=[-500 -200];
[~,t_start_b]=min(abs(scl.t_ms-baseline_region(1)));
[~,t_end_b]=min(abs(scl.t_ms-baseline_region(2)));

event_region=[200 400];
[~,t_start]=min(abs(scl.t_ms-event_region(1)));
[~,t_end]=min(abs(scl.t_ms-event_region(2)));

axis_lims=[-4 3];
spec_chans=[7 25 58];
unit_string=' ln(uV^2)';

for hyp=1:3
figure;
pair_hypinds = opt.pair_inds == hyp;
pair_hypchans=unique(opt.coherence_pairs(pair_hypinds,:));
%pair_hypchans=spec_chans(hyp);
overtitle=opt.pair_indlbls{hyp};
for freq_range=1:4
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    %[~,f_ind]=min(abs(scl.freqs-pp.f_indiv_hz(freq_range)));
    subplot(1,4,freq_range);
    subtitle = sprintf('%1.1f - %1.1f Hz',scl.freqs(f_start), ...
        scl.freqs(f_end) );
    hold on;
for group=pp.chosen_g(pp.plotn_g)
    
    x_scatterdata = squeeze(meanx( log( rpowdata(pair_hypchans, f_end:f_start, s_inds_g_r(:,group))),3));
    y_scatterdata = meanx( log( wave_totpowdata(t_start_b:t_end_b, pair_hypchans, f_end:f_start, :, s_inds_g(:,group))), 5);
    z_scatterdata = meanx( log( wave_totpowdata(t_start:t_end, pair_hypchans, f_end:f_start, :, s_inds_g(:,group))), 5);
    
    %x_scatterdata = squeeze(meanx(rcohdata(pair_hypinds, f_ind, s_inds_g_r(:,group)),3));
    %y_scatterdata = meanx( cohdata(t_start_b:t_end_b, f_ind, :, pair_hypinds, s_inds_g(:,group)), 5);
    
    scatter3( x_scatterdata, y_scatterdata, z_scatterdata, scl.g_color{group} );
    
    [x,y,z] = fitline_3dpts(x_scatterdata, y_scatterdata, z_scatterdata, axis_lims.*2);
    plot3(x,y,z,scl.g_color{group});
    
end
axis([axis_lims(1) axis_lims(2) axis_lims(1) axis_lims(2) axis_lims(1) axis_lims(2)]); grid on;
xlabel(['Resting',unit_string]); ylabel(['Baseline',unit_string]); zlabel(['Event',unit_string]);
title(subtitle);
plot3(linspace(axis_lims(1),axis_lims(2),100),linspace(axis_lims(1),axis_lims(2),100), ...
    linspace(axis_lims(1),axis_lims(2),100),'k--');
hold off;
%view(-20, 25);
%view(90, 0); %baseline = x, event = y
view(0, 90); %resting = x, baseline = y
end
%tightfig;
plottitle(overtitle,2);
set(gcf, 'Position', [1300 100 1900 325]);
end