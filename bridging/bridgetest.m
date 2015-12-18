function bridgetest(in_list,out_log,condset,subset)
% given a list of *.avg.h1 filepaths, generates output diagnostic for the
% presence of electrode bridging or other issues in the recording.
% outputs are 'in_list-log.txt' and 'in_list.pdf' which will appear in the
% current directory

% example inputs:
% in_list='/export/home/mike/matlab/masterlist.txt';
% out_log='/export/home/mike/matlab/masterlist-log.txt';
% condset=1;
% subset=[6 14 22];


h1_list = list2cell(in_list); % import the list as a cell array

fid = fopen(out_log, 'w'); %open the output log for writing
fprintf(fid,'Bridging test run on %s\n',datestr(clock));
fprintf(fid,'---\n\n');

tempdir = [cd,filesep,'bridgetemp'];
if ~exist(tempdir,'dir')
    mkdir(tempdir)
end

if nargin < 4
    subset = 1:length(h1_list);
end

alpha = .001; %these are default for now
chanlocs_file = '/active_projects/matlab_common/61chans_ns.mat';
load(chanlocs_file); %as 'chan_locs'

figure('Visible','Off');
pdfdummy=0;

for rec=subset %for each separate *.avg.h1 file

[~,name] = fileparts(h1_list{rec});
fprintf(fid,'Subject/session:\t\t%s\n',name);
    
clf;
    
h1_struct = read_hdf1_dataV7(h1_list{rec});
h1_data = h1_struct.data_struct.hdf1_avg_data(:,[1:31,33:62],:);
%chan_labels = h1_struct.run_struct.channel_label([1:31,33:62]);
[n_conds,n_chans,~] = size(h1_data);

if nargin < 3
    condset = 1:n_conds;
end

for cond=condset %for each condition

cond_name = h1_struct.case_struct.descriptor{cond};
fprintf(fid,'Condition:\t\t\t\t%s\n',cond_name);
    
erpdata = squeeze(h1_data(cond,:,:))'; %erpdata is timepts x channels

% correlation matrix
rho = corr(erpdata);
subplot(2,2,1);
imagesc(rho,[0.4 1]);
colorbar;
freezeColors;

% determining mean correlation of each channel and outliers
meancorr = mean(rho - eye(size(rho)));
[~,outlier_inds] = deleteoutliers(meancorr,alpha);
gmeancorr = mean(meancorr);

% histogram of channel correlations
subplot(2,2,2);
histx(squareform(tril(rho,-1)));
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .8 .8]);

% topoplot of mean correlations
subplot(2,2,3);
maplims = [min(meancorr) max(meancorr)];
topoplot(meancorr,chan_locs,'maplimits',maplims,'colormap',jet);
colorbar;

% PCA of data - bridged data will require fewer comps to account for 90%
subplot(2,2,4);
s = svd(erpdata);
svd_curve = cumsum(s)/sum(s);
plot(svd_curve);
axis([1 n_chans 0 1]);
grid on;
crit = find(svd_curve>.9,1);
text(20,0.6,sprintf('MeanCorr = %0.2f',gmeancorr));
text(20,0.4,sprintf('CritComps = %d',crit));

% title for figure
overtitle = sprintf('%s - %s',name,cond_name);
plottitle(overtitle,2);

% save the figure as a pdf
figname = [name,'_',cond_name,'.pdf'];
pdfdummy=pdfdummy+1;
pdfpath{pdfdummy} = fullfile(tempdir,figname);
print(pdfpath{pdfdummy},'-dpdf')

% text ouput
fprintf(fid,'Mean correlation:\t\t%0.2f\n',gmeancorr);
fprintf(fid,'# of comps to 90%%:\t\t%d\n',crit);
fprintf(fid,'Potentially bad chans:\t');
if ~isempty(outlier_inds)
    for chan=outlier_inds
        fprintf(fid,'%s, ',chan_locs(chan).labels);
    end
else
    fprintf(fid,'none');
end
fprintf(fid,'\n\n');

end %condition loop
end %subject loop

% append all of the pdfs
[~,list_name]=fileparts(in_list);
append_pdfs([list_name,'.pdf'],pdfpath{:});

% delete the temporary stuff
delete([tempdir,filesep,'*.pdf']);

if fid > 1
    fclose(fid);
end

end %function