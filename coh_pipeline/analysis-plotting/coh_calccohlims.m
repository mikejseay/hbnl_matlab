%% calculation part 1 - making matrices of min/max, mean, and stds

s_inds=find(s_inds_g(:,9))';
cohlim_mat=zeros(2,length(s_inds), length(pp.plotn_f), imp.maxpairs);
cohmean_mat=zeros(length(s_inds), length(pp.plotn_f), imp.maxpairs);
cohstd_mat=zeros(length(s_inds), length(pp.plotn_f), imp.maxpairs);
sdum=0;

timeregion=[-200 800];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));

for s=s_inds
sdum=sdum+1;
for freq_range=pp.plotn_f
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
for pair=1:imp.maxpairs
    cohlim_mat(:,sdum,freq_range,pair) = ...
        [minx(cohdata(t_start:t_end,f_end:f_start,:,pair,s),[]) ...
        maxx(cohdata(t_start:t_end,f_end:f_start,:,pair,s),[]) ];
    cohmean_mat(sdum,freq_range,pair) = ...
        meanx(cohdata(t_start:t_end,f_end:f_start,:,pair,s),[]);
    cohstd_mat(sdum,freq_range,pair) = ... 
        std(meanx(cohdata(t_start:t_end,f_end:f_start,:,pair,s),1));
    %here we mean std across time - how much does the coherence express
    % variation throughout the event range
end
end
end

%% calculation but better

s_inds=find(s_inds_g(:,9))';
p_dists=pair_distance(opt.coherence_pairs,chan_locs);

cohmean_mat=permute(meanx(cohdata,[2 4 5]),[3 1 2]);
%cohmean_mat=zeros(
%for freq_range=1:pp.plotn_f
%    cohmean_mat(:,freq_range,:)=permute(meanx(cohdata,[2 4 5]),[3 1 2]);
%end

cohmin_mat=permute(minx(cohdata,[2 4 5]),[3 1 2]);
cohmax_mat=permute(maxx(cohdata,[2 4 5]),[3 1 2]);
cohlim_mat(1,:,:,:)=cohmin_mat;
cohlim_mat(2,:,:,:)=cohmax_mat;

%% plotting part 1 - scatter mean coherence vs distance for each range

figure;
for freq_range=pp.plotn_f
    subplot(2,4,freq_range)
    
    scatter(p_dists,meanx(cohmean_mat(:,freq_range,:),3))
    axis([min(p_dists) max(p_dists) 0 1]);
    
    title(sprintf('%1.1f - %1.1f Hz',pp.f_start_hz(freq_range),pp.f_end_hz(freq_range)))
    
end

%% plotting part 1.2 - boxplot coherence for each pair

[p_dists_sorted,inds]=sort(p_dists);

pair_types={'intra-frontal','intra-parietal','intra-occipital','frontal-parietal',...
    'frontal-occipital','occipito-parietal'};
pair_type_inds={1:13,14:26,27:37,38:59,60:77,78:92};

pair_g=zeros(size(p_dists));
for ptype=1:length(pair_types)
    pair_g(pair_type_inds{ptype})=ptype;
end

for freq_range=3
    figure;
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    
    boxplot(meanx(cohmean_mat(:,f_end:f_start,inds),[1 3]),'notch','on',...
        'labels',scl.p_label,'labelorientation','inline');
    axis([0 length(p_dists)+1 0 1]);
    %boxplot(meanx(cohmean_mat(:,freq_range,:),[1 3]),pair_g,'notch','on',...
    %    'labels',pair_types)
    %axis([min(p_dists) max(p_dists) 0 1]);
    
    title(sprintf('%1.1f - %1.1f Hz',pp.f_start_hz(freq_range),pp.f_end_hz(freq_range)))
    
end

% use 'colorgroup' for the sets of pairs (frontal/parietal/occipital)
% 

%% plotting part 1.3 - boxplots of condition differences over time



%% calc pt 2 - negative exponential fit



%% plotting part 2 - scatter coherence for each subject separately

freq_range=3;

figure;
sdum=0;
for s=s_inds
    sdum=sdum+1;
    subplot(6,10,sdum)
    
    scatter(p_dists,squeeze(cohmean_mat(sdum,freq_range,:)));
    
end

%% plotting part 3 - scatter coherence +/- SD vs distance

figure;
for freq_range=pp.plotn_f
    subplot(2,4,freq_range)
    coh_plotmeans=meanx(cohmean_mat(:,freq_range,:),3);
    %scatter(p_dists,coh_plotmeans,'.'); hold on;
    for pair=1:imp.maxpairs
        coh_plotstd=meanx(cohstd_mat(:,freq_range,pair),[]);
        line([p_dists(pair) p_dists(pair)], ...
            [coh_plotmeans(pair) - coh_plotstd, coh_plotmeans(pair) + coh_plotstd], ...
            'Marker','+','Color',scl.p_color(pair,:)); hold on;
    end
    hold off;
    
    title(sprintf('%1.1f - %1.1f',pp.f_start_hz(freq_range),pp.f_end_hz(freq_range)))
    
end
