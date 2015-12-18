function bridgetest_old(in_list,cond,subset)
% given a list of *.avg.h1 filepaths, generates output diagnostic for the
% presence of electrode bridging or other issues in the recording.

% example inputs:
%in_list='/export/home/mike/matlab/masterlist.txt';
%subset=[6 14 22];
%cond=1;

h1_list=list2cell(in_list);

if nargin<3
    subset=1:length(h1_list);
end
if nargin<2
    cond=1;
end

spacer=3;
alpha=.001;
chanlocs_file='/active_projects/matlab_common/61chans_ns.mat';

load(chanlocs_file); %as 'chan_locs'

figure;

for rec=subset

clf;
    
h1_struct=read_hdf1_dataV7(h1_list{rec});

h1_data=h1_struct.data_struct.hdf1_avg_data(:,[1:31,33:62],:);
chan_labels=h1_struct.run_struct.channel_label([1:31,33:62]);

[n_conds,n_chans,n_timepts]=size(h1_data);

x = squeeze(h1_data(cond,:,:))';

rho = corr(x);

s1=subplot(221);
imagesc(rho,[0.4 1]);
set(gca,'XTick',1:spacer:n_chans,'YTick',1:spacer:n_chans,'XTickLabel',chan_labels(1:spacer:end), ...
    'YTickLabel',chan_labels(1:spacer:end));
freezeColors;

meancorr=mean(rho-eye(size(rho)));
[~,outlier_inds]=deleteoutliers(meancorr,alpha);
gmeancorr=mean(meancorr);

subplot(222);
%hist(squareform(tril(rho,-1)),50);
histx(squareform(tril(rho,-1)));
h = findobj(gca,'Type','patch');
h.FaceColor = [.8 .8 .8];

subplot(223);
%bar(meancorr, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor', [1 1 1]);
%axis([1 n_chans 0 1]);
%set(gca,'XTick',1:spacer:n_chans,'XTickLabel',chan_labels(1:spacer:end));
maplims=[min(meancorr) max(meancorr)];
topoplot(meancorr,chan_locs,'maplimits',maplims,'colormap',parula);
%colorbar;

subplot(224);
s=svd(x);
svd_curve=cumsum(s)/sum(s);
plot(svd_curve);
axis([1 n_chans 0 1]);
grid on;
crit=find(svd_curve>.9,1);
text(20,0.6,sprintf('MeanCorr = %0.2f',gmeancorr));
text(20,0.4,sprintf('CritComps = %d',crit));

[~,name]=fileparts(h1_list{rec});
cond_name=h1_struct.case_struct.descriptor{cond};
overtitle=sprintf('%s - %s',name,cond_name);
plottitle(overtitle);

%fprintf('Grand mean correlation was %0.2f\n',gmeancorr);
%fprintf('Number of components accounting for 90%% variance was %d\n',crit);
fprintf('Potentially abnormal channels: ');
if ~isempty(outlier_inds)
    for chan=outlier_inds
        fprintf('%s, ',chan_locs(chan).labels);
    end
else
    fprintf('none')
end
fprintf('\n');

end

end