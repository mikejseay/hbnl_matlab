%% calculating baseline exponential fit for all subjects/frequency ranges

s_inds=find(s_inds_g(:,1))';

p_dists=pair_distance(opt.coherence_pairs,chan_locs);
a=min(p_dists); b=max(p_dists);
distance=linspace(a,b,100);

timeregion=[-500 -200];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));

n_params=3;
params_mat=zeros(n_params,length(pp.plotn_f),length(s_inds));
error_mat=zeros(length(pp.plotn_f),length(s_inds));
estimates_mat=zeros(length(p_dists),length(pp.plotn_f),length(s_inds));

%options = optimset('MaxFunEvals',1000);
options = optimset('MaxFunEvals',[]);

for s=s_inds
for freq_range=pp.plotn_f

    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));

    fitdata=meanx(cohdata(t_start:t_end,f_end:f_start,:,:,s),4);

    [params_mat(:,freq_range,s),model]=exp_fit3(fitdata,p_dists,options);
    [error_mat(freq_range,s),estimates_mat(:,freq_range,s)]= ...
        model(params_mat(:,freq_range,s));

end
end

%% for a given subject/frequency plot the fit

timeregion=[-500 -200];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));

g_lab=char(s_demogs.group(s));
subplot_dummy=0;
for s=s_inds
for freq_range=3 %pp.plotn_f
    if mod(subplot_dummy,15)==0
        figure;
    end
    subplot_dummy=subplot_dummy+1;
    subplot(3,5,mod(subplot_dummy,15)+1);

    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));

    plotdata=meanx(cohdata(t_start:t_end,f_end:f_start,:,:,s),4);

    model_line=exp(params_mat(1,freq_range,s)+params_mat(2,freq_range,s)* ...
        distance)+params_mat(3,freq_range,s);
    
    scatter(p_dists,plotdata); hold on;
    plot(distance, model_line); hold off;
    
    text(.75*max(p_dists), .75, sprintf('%1.1f',error_mat(freq_range,s)));
    
    title(sprintf('S %d / %s',s,g_lab(1:3)));
end
end

%% for a given subject/frequency, plot the event-related condition means relative to baseline fit

sstyles={'gs','rs'};

timeregion=[200 500];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));

subplot_dummy=0;
for s=s_inds
for freq_range=1:2 %pp.plotn_f
    if mod(subplot_dummy,15)==0
        overtitle=sprintf('%d - %d ms / %1.1f - %1.1f Hz',timeregion(1), ...
            timeregion(2), pp.f_start_hz(freq_range), pp.f_end_hz(freq_range));
        plottitle(overtitle);
        axis_prunelabels('xy');
        tightfig;
        set_print_size(20,20/scl.phi);
        figure;
    end
    subplot_dummy=subplot_dummy+1;
    subplot(3,5,mod(subplot_dummy,15)+1);

    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
    
    for cond=1:2
    plotdata(:,cond)=meanx(cohdata(t_start:t_end,f_end:f_start,cond,:,s),4);
    scatter(p_dists,plotdata(:,cond),sstyles{cond},'SizeData',10); hold on;
    end
    for pair=1:imp.maxpairs
    line([p_dists(pair) p_dists(pair)], [plotdata(pair,1) plotdata(pair,2)], ...
        'Color','k'); hold on;
    end
    
    model_line=exp(params_mat(1,freq_range,s)+params_mat(2,freq_range,s)* ...
        distance)+params_mat(3,freq_range,s);
    plot(distance, model_line); hold off;
    axis([a-1 b+1 0 1]); grid on;
    g_lab=char(s_demogs.group(s));
    title(sprintf('S %d / %s',s,g_lab(1:3)));
end
end
plottitle(overtitle);
axis_prunelabels('xy');
tightfig;
set_print_size(20,20/scl.phi);

%% subject-mean fit with data cursor labels for pairs

freqregion=[4 7];
[~,f_start]=min(abs(scl.freqs-freqregion(1)));
[~,f_end]=min(abs(scl.freqs-freqregion(2)));


%fit
timeregion=[-500 -200];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));

fitdata=meanx(cohdata(t_start:t_end,f_end:f_start,:,:,s_inds_g(:,1)),4);
%fitdata=median(meanx(cohdata(t_start:t_end,f_end:f_start,:,:,s_inds_g(:,1)),[4 5]),2);
[params,model]=exp_fit3(fitdata,p_dists,options);
[error,estimates]=model(params);
model_line=exp(params(1)+params(2)*distance)+params(3);

%plot
%timeregion=[200 500];
timeregion=[-500 -200];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));
plotdata=meanx(cohdata(t_start:t_end,f_end:f_start,:,:,s_inds_g(1,:)),4);

global lookup scl

lookup=[p_dists,plotdata];

f=figure;
plot(distance, model_line); hold on;
scatter(p_dists, plotdata); hold off;
dcm_h=datacursormode(f);
set(dcm_h,'UpdateFcn',@coh_pairlabel_dcm)

%% subject-mean fit with separate colors for each hypothesis

timeregion=[-500 -200];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));

hyps=length(opt.pair_indlbls);
hyp_colors=distinguishable_colors(hyps);

lookup=[];

for freq_range=2 %pp.plotn_f
    
    [~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
    [~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));

    plotdata=meanx(cohdata(t_start:t_end,f_end:f_start,:,:,s_inds_g(1,:)),4);
    
    [params,model]=exp_fit3(plotdata,p_dists,options);
    [error,estimates]=model(params);
    model_line=exp(params(1)+params(2)*distance)+params(3);

    lookup=[p_dists,plotdata];
    
    f=figure;
    overtitle=sprintf('%1.1f - %1.1f Hz, %d - %d ms',pp.f_start_hz(freq_range), ...
        pp.f_end_hz(freq_range), timeregion(1), timeregion(2));
    
    for hyp=1:hyps
        plot_hypinds=opt.pair_inds==hyp;
        scatter(p_dists(plot_hypinds), plotdata(plot_hypinds),'MarkerEdgeColor',hyp_colors(hyp,:));
        hold on;
    end
    legend(opt.pair_indlbls)
    plot(distance, model_line); hold off;
    axis([a-1 b+1 0 1]);
    dcm_h=datacursormode(f);
    set(dcm_h,'UpdateFcn',@coh_pairlabel_dcm)
    
    xlabel('Distance (cm)');
    ylabel('ERPCOH');
    
    plottitle(overtitle);
    
end

%% difference from fit in same way as above

lookup=[p_dists,plotdata - estimates];

f=figure;
scatter(p_dists, plotdata - estimates); hold on;
hline(0,'k--'); hold off;
dcm_h=datacursormode(f);
set(dcm_h,'UpdateFcn',@coh_pairlabel_dcm)

% why don't we see negative difference from fit for long pairs???
% maybe it's because the fit is better for long-range pairs

%% percent difference (division) from fit in same way as above

lookup=[p_dists,plotdata ./ estimates];

f=figure;
scatter(p_dists, plotdata ./ estimates); hold on;
hline(1,'k--'); hold off;
dcm_h=datacursormode(f);
set(dcm_h,'UpdateFcn',@coh_pairlabel_dcm)

%% condition difference from in same way as above (SHOULD BE SAME AS 2 BELOW)

%basetimeregion=[

timeregion=[200 400];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));
freqregion=[4 7];
[~,f_start]=min(abs(scl.freqs-freqregion(1)));
[~,f_end]=min(abs(scl.freqs-freqregion(2)));

lookup=[];

f=figure;
overtitle=sprintf('%1.1f - %1.1f Hz',freqregion(1),freqregion(2));
subplot_dummy=0;
for group=pp.chosen_g(pp.plotn_g)
%subplot_dummy=subplot_dummy+1;
%subplot(1,2,subplot_dummy);
plotdata=zeros(length(p_dists),2);
sstyles={'gs','rs'};
for cond=1:2
    plotdata(:,cond)=meanx(cohdata(t_start:t_end,f_end:f_start,cond,:,s_inds_g(:,group)),4);
    lookup=[lookup;p_dists,plotdata(:,cond) - estimates];
    scatter(p_dists, plotdata(:,cond) - estimates,sstyles{cond}); hold on;
end
for pair=1:imp.maxpairs
    line([p_dists(pair) p_dists(pair)], [plotdata(pair,1) - estimates(pair) ...
        plotdata(pair,2) - estimates(pair)],'Color','k'); hold on;
end
hline(0,'k--'); hold off;
axis([a-1 b+1 -0.15 0.4]);
xlabel('Distance (cm)'); ylabel('ERPCOH minus baseline fit');
dcm_h=datacursormode(f);
set(dcm_h,'UpdateFcn',@coh_pairlabel_dcm)
title(scl.g_label{group})
end
plottitle(overtitle);
axis_prunelabels('xy');
tightfig;
set_print_size(20,20/scl.phi);

%% condition difference from in same way as above but per subject/freq estimates

%basetimeregion=[

timeregion=[200 500];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));
freqregion=[4 7];
[~,f_start]=min(abs(scl.freqs-freqregion(1)));
[~,f_end]=min(abs(scl.freqs-freqregion(2)));

%lookup=[];

subplot_dummy=0;
for freq_range=1:2
for s=s_inds
if mod(subplot_dummy,15)==0
    tightfig;
    figure;
end
subplot_dummy=subplot_dummy+1;
subplot(3,5,mod(subplot_dummy,15)+1);
%overtitle=sprintf('S %d / %1.1f - %1.1f Hz',s,freqregion(1),freqregion(2));

plotdata=zeros(length(p_dists),2);
sstyles={'gs','rs'};
for cond=1:2
    plotdata(:,cond)=meanx(cohdata(t_start:t_end,f_end:f_start,cond,:,s),4);
    %lookup=[lookup;p_dists,plotdata(:,cond) - estimates_mat(:,freq_range,s)];
    scatter(p_dists, plotdata(:,cond) - estimates_mat(:,freq_range,s),sstyles{cond}); hold on;
end
for pair=1:imp.maxpairs
    line([p_dists(pair) p_dists(pair)], [plotdata(pair,1) - estimates_mat(pair,freq_range,s) ...
        plotdata(pair,2) - estimates_mat(pair,freq_range,s)],'Color','k'); hold on;
end
hline(0,'k--'); hold off;
axis([a-1 b+1 -0.5 0.7]);
%xlabel('Distance (cm)'); ylabel('ERPCOH minus baseline fit');
%dcm_h=datacursormode(f);
%set(dcm_h,'UpdateFcn',@coh_pairlabel_dcm)
%plottitle(overtitle);
end
end
tightfig;

%% condition difference in same way as above but use per subject/freq estimates
%  to normalize each before taking the mean...
% THIS SHOULD BE EXACT SAME AS 2 ABOVE

timeregion=[200 500];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));

for freq_range=1:2
f=figure;
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
overtitle=sprintf('%1.1f - %1.1f Hz',pp.f_start_hz(freq_range), ...
    pp.f_end_hz(freq_range));
%subplot_dummy=0;
for group = pp.chosen_g(pp.plotn_g)
%subplot_dummy=subplot_dummy+1;
%subplot(1,2,subplot_dummy);
plotdata=zeros(length(p_dists),2);
sstyles={'gs','rs'};
for cond=1:2
    plotdata(:,cond)=mean( meanx(cohdata(t_start:t_end,f_end:f_start,cond,:,s_inds_g(:,group)),[4 5]) - ...
        squeeze(estimates_mat(:,freq_range,s_inds_g(:,group))) ,2);
    %plotdata(:,cond)=mean( meanx(cohdata(t_start:t_end,f_end:f_start,cond,:,s_inds_g(:,group)),[4 5]) ./ ...
    %    squeeze(estimates_mat(:,freq_range,s_inds_g(:,group))) ,2);
    %lookup=[lookup;p_dists,plotdata(:,cond) - estimates_mat(:,freq_range,s)];
    scatter(p_dists, plotdata(:,cond),sstyles{cond}); hold on;
end
for pair=1:imp.maxpairs
    line([p_dists(pair) p_dists(pair)], [plotdata(pair,1) plotdata(pair,2)], ...
        'Color','k'); hold on;
end
hline(0,'k--'); hold off;
axis([a-1 b+1 -0.1 0.4]);
xlabel('Distance (cm)'); ylabel('ERPCOH minus baseline fit');
title(scl.g_label{group})
%dcm_h=datacursormode(f);
%set(dcm_h,'UpdateFcn',@coh_pairlabel_dcm)
end
plottitle(overtitle);
tightfig;
end

%% condition difference in same way as above but use per subject/freq estimates
%  to divisively normalize each before taking the mean...
% THIS IS NOT THE EXACT SAME AS 2 ABOVE

timeregion=[200 400];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));

hyp_colors=distinguishable_colors(length(opt.pair_indlbls),'g');
clear overtitle;
for freq_range=1:3
f=figure;
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
overtitle=sprintf('%1.1f - %1.1f Hz, %d - %d ms',pp.f_start_hz(freq_range), ...
    pp.f_end_hz(freq_range), timeregion(1), timeregion(2));
%subplot_dummy=0;
for group = pp.chosen_g(pp.plotn_g)
%subplot_dummy=subplot_dummy+1;
%subplot(1,2,subplot_dummy);
plotdata=zeros(length(p_dists),2);
sstyles={'gs','rs'};
for cond=1:2
    plotdata(:,cond)=mean( meanx(cohdata(t_start:t_end,f_end:f_start,cond,:,s_inds_g(:,group)),[4 5]) ./ ...
        squeeze(estimates_mat(:,freq_range,s_inds_g(:,group))) ,2);
    lookup=[lookup;p_dists,plotdata(:,cond)];
    scatter(p_dists, plotdata(:,cond),sstyles{cond}); hold on;
end
for pair=1:imp.maxpairs
    line([p_dists(pair) p_dists(pair)], [plotdata(pair,1) plotdata(pair,2)], ...
        'Color',hyp_colors(opt.pair_inds(pair),:),'LineWidth',2); hold on;
end
hline(1,'k--'); hold off;
axis([a-1 b+1 0.7 2.1]);
xlabel('Distance (cm)'); ylabel('ERPCOH divided by baseline fit');
title(scl.g_label{group});

for hyp=1:length(opt.pair_indlbls)
    text((.75)*max(p_dists), 2-(hyp/20), sprintf('%s', opt.pair_indlbls{hyp} ),'Color', ...
        hyp_colors(hyp,:)); hold on;
end

dcm_h=datacursormode(f);
set(dcm_h,'UpdateFcn',@coh_pairlabel_dcm)
end
plottitle(overtitle);
tightfig;
end

%% condition difference in same way as above but use per subject/freq estimates
%  to divisively normalize each to BASELINE COHERENCE before taking the mean...

baseline_region=[-500 -200];
[~,t_start_b]=min(abs(scl.t_ms-baseline_region(1)));
[~,t_end_b]=min(abs(scl.t_ms-baseline_region(2)));

timeregion=[200 400];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));

hyp_colors=distinguishable_colors(length(opt.pair_indlbls),'g');
clear overtitle;
for freq_range=1:3
f=figure;
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));
overtitle=sprintf('%1.1f - %1.1f Hz, %d - %d ms',pp.f_start_hz(freq_range), ...
    pp.f_end_hz(freq_range), timeregion(1), timeregion(2));
%subplot_dummy=0;
for group = pp.chosen_g(pp.plotn_g)
%subplot_dummy=subplot_dummy+1;
%subplot(1,2,subplot_dummy);
plotdata=zeros(length(p_dists),2);
sstyles={'gs','rs'};
for cond=1:2
    plotdata(:,cond)=mean( meanx(cohdata(t_start:t_end,f_end:f_start,cond,:,s_inds_g(:,group)),[4 5]) ./ ...
        meanx(cohdata(t_start_b:t_end_b,f_end:f_start,:,:,s_inds_g(:,group)),[4 5]) ,2);
    lookup=[lookup;p_dists,plotdata(:,cond)];
    scatter(p_dists, plotdata(:,cond),sstyles{cond}); hold on;
end
for pair=1:imp.maxpairs
    line([p_dists(pair) p_dists(pair)], [plotdata(pair,1) plotdata(pair,2)], ...
        'Color',hyp_colors(opt.pair_inds(pair),:),'LineWidth',2); hold on;
end
hline(1,'k--'); hold off;
axis([a-1 b+1 0.7 2.1]);
xlabel('Distance (cm)'); ylabel('ERPCOH divided by baseline');
title(scl.g_label{group});

for hyp=1:length(opt.pair_indlbls)
    text((.75)*max(p_dists), 2-(hyp/20), sprintf('%s', opt.pair_indlbls{hyp} ),'Color', ...
        hyp_colors(hyp,:)); hold on;
end

dcm_h=datacursormode(f);
set(dcm_h,'UpdateFcn',@coh_pairlabel_dcm)
end
plottitle(overtitle);
tightfig;
end


%% legacy calculation but better

s_inds=find(s_inds_g(:,1))';
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

chan_cohlims=meanx(cohlim_mat,[1 4])';

%% plotting each range by distance
p_dists=pair_distance(opt.coherence_pairs,chan_locs);
figure;
for pair=1:size(chan_cohlims,1)
    line([p_dists(pair) p_dists(pair)],chan_cohlims(pair,:),'Marker','+','Color',scl.p_color(pair,:)); hold on;
end
ylabel('Coherence Min/Max')
xlabel('Distance (cm)')

%% scattering pair distance against coherence range
chan_cohranges=diff(chan_cohlims,1,2);
figure; scatter(p_dists,chan_cohranges);
xlabel('Distance (cm)')
ylabel('Coherence Range')

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

for freq_range=6
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

%% plotting part 1.3 - boxplots of condition differences in time windows

timeregion=[200 500];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));
freqregion=[4 5.3];
[~,f_start]=min(abs(scl.freqs-freqregion(1)));
[~,f_end]=min(abs(scl.freqs-freqregion(2)));

[p_dists_sorted,inds]=sort(p_dists);

pair_types={'intra-frontal','intra-parietal','intra-occipital','frontal-parietal',...
    'frontal-occipital','occipito-parietal'};
pair_type_inds={1:13,14:26,27:37,38:59,60:77,78:92};

pair_g=zeros(size(p_dists));
for ptype=1:length(pair_types)
    pair_g(pair_type_inds{ptype})=ptype;
end

figure;
plotdata = meanx(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{1},:,s_inds_g(:,1)),[4 5]) - ...
    meanx(cohdata(t_start:t_end,f_end:f_start,pp.cond_diff{2},:,s_inds_g(:,1)),[4 5]);
boxplot(plotdata(inds,:)','notch','on',...
    'labels',{scl.p_label{inds}},'labelorientation','inline');
axis([0 length(p_dists)+1 -.5 .5]);
hold on; hline(0,'k--'); hold off;

%%

% bull shit follows
%{ 
beta=regress(log(data),[ones(length(p_dists),1) p_dists]);

figure;
scatter(p_dists,log(data)); hold on;
plot(linspace(a,b,100),linspace(a,b,100)*beta(2)+beta(1));



tp=[0.8, 1, -2, 0.2];
my_f=@(x)tp(1).^(tp(2)*(x+tp(3)))+tp(4);
xvals=linspace(a,b,100);
yvals=my_f(xvals);
figure; plot(xvals,yvals);


tp=[-0.5 -2 0.2];
my_f=@(x)exp(tp(1)*(x+tp(2)))+tp(3);
xvals=linspace(a,b,100);
yvals=my_f(xvals);
%figure;
hold on; plot(xvals,yvals); hold off;
axis([a b 0 1]);



[p_dists_s,inds]=sort(p_dists,'descend');
figure; scatter(p_dists_s,data(inds));


[k, yInf, y0, yFit] = fitExponential(p_dists_s, data(inds));

xvals = linspace(a,b,100);
yvals = yInf + (y0-yInf) * exp(-k*(xvals-a));
hold on; plot(xvals,yvals); hold off;

%% plot mean coherence vs. mean baseline fit in TF windows for each condition

sp_rowlabel=make_freqlabels(pp.f_start_hz(pp.plotn_f),pp.f_end_hz(pp.plotn_f));
sp_columnlabel=make_timelabels(pp.t_start_ms(pp.plotn_t),pp.t_end_ms(pp.plotn_t));
%sp_columnlabel{end}='';
x_plotlabel='Time Windows';
y_plotlabel='Frequencies';

timeregion=[-500 -200];
[~,t_start]=min(abs(scl.t_ms-timeregion(1)));
[~,t_end]=min(abs(scl.t_ms-timeregion(2)));
fitdata=meanx(cohdata(t_start:t_end,:,:,:,s_inds_g(1,:)),4);
[params,~]=exp_fit3(fitdata,p_dists);
model_line=exp(params(1)+params(2)*distance)+params(3);

for cond=1:2
figure;
subplot_dum=0;
overtitle=scl.cond_label{cond};
for freq_range=pp.plotn_f
[~,f_start]=min(abs(scl.freqs-pp.f_start_hz(freq_range)));
[~,f_end]=min(abs(scl.freqs-pp.f_end_hz(freq_range)));

for win=pp.plotn_t
[~,t_start]=min(abs(scl.t_ms-pp.t_start_ms(win)));
[~,t_end]=min(abs(scl.t_ms-pp.t_end_ms(win)));

subplot_dum=subplot_dum+1;
subplot(length(pp.plotn_f),length(pp.plotn_t),subplot_dum)

data=meanx(cohdata(t_start:t_end,f_end:f_start,:,:,s_inds_g(:,1)),4);

scatter(p_dists,data,'.'); hold on; %actual data

plot(distance,model_line); hold off; %full model line
axis([a-1 b+1 0 1]);

end
end
adorn_plots(sp_rowlabel,sp_columnlabel,x_plotlabel,y_plotlabel,overtitle,[length(pp.plotn_f),length(pp.t_start_ms)]);
end

%}