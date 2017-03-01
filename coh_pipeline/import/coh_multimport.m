function [matlist_all, imp, n_trials_all, behdata, erpdata, wave_totpowdata, ...
    wave_evknormdata, cohdata, phidata, wave_totdata] = ...
    coh_multimport(opt, demogsfile, tf_downsamp, behmat_path, datatypes) %wave_evknorm
% Imports computed EEG measures into the MATLAB workspace for freestyle
% plotting and exploration.

n_opt = length(opt);

if nargin < 3
    tf_downsamp=2; % downsampling by half is default
end

% make a list of each subject's .mat file
if nargin < 2 % if no demogsfile was specified, use this
    matlist_all=coh_updatemats(opt(1).outpath);
else
    matlist_all = cell(1, n_opt);
    for o = 1:n_opt
        matlist_all{o} = coh_updatemats(opt(o).outpath, demogsfile);
    end
    % attempt a lateral merge
    matlist_all = merge_cellcell( matlist_all );
end

% if a behavioral path has been specified, and behavioral data has been
% requested
if nargin>=4 && nargout >=4 && ~isempty(behmat_path)
    behdata=beh_import(behmat_path,matlist_all);
else
    behdata=[];
end

ns=size(matlist_all,1);

%size of data dimensions, taking into account a time-downsampling factor
erpmaxtimepts = opt(1).n_samps - 4;
times2import=1:erpmaxtimepts;
erptimerate = opt(1).rate;

if isfield(opt, 'tf_timedownsamp_ratio')
    tf_downsamp = tf_downsamp / opt(1).tf_timedownsamp_ratio;
    tfmaxtimepts = round((opt(1).n_samps)/(opt(1).tf_timedownsamp_ratio*tf_downsamp)) - 1;
    tf_endpt = round(opt(1).n_samps/opt(1).tf_timedownsamp_ratio) - 1;
    tftimes2import=1:tf_downsamp:tf_endpt;
    tftimerate=opt(1).rate/(tf_downsamp*opt(1).tf_timedownsamp_ratio);
else
    tfmaxtimepts = round((opt(1).n_samps)/tf_downsamp);
    tf_endpt = opt(1).n_samps - 1;
    tftimes2import=1:tf_downsamp:tf_endpt;
    tftimerate=opt(1).rate/tf_downsamp;
end

maxfreqs=length(opt(1).wavelet_scales);
%maxconds=length(opt(1).case_vec);
maxconds = sum(cellfun(@length, {opt.case_vec}));
maxchans=length(opt(1).chan_vec);
maxpairs=size(opt(1).coherence_pairs,1);

%indicate mats to load
if nargin<5
    n_datatypes=min(nargout-2, 5); %6
    datatype_names={'n_trials','erp', 'wave_totpow', 'wave_evknorm', 'coh'};
    datatypes={datatype_names{1:n_datatypes}};
    %if nargout>=10
    %    datatypes{end+1}='wave_evknorm';
    %end
end

%pre-allocate vars
n_trials_all=zeros(maxconds,ns);
erpdata=zeros(erpmaxtimepts,maxchans,maxconds,ns);
%wave_evkdata=zeros(maxtimepts,maxchans,maxfreqs,maxconds,ns);
wave_evknormdata=zeros(tfmaxtimepts,maxchans,maxfreqs,maxconds,ns);
wave_totdata=zeros(tfmaxtimepts,maxchans,maxfreqs,maxconds,ns);
wave_totpowdata=zeros(tfmaxtimepts,maxchans,maxfreqs,maxconds,ns);
cohdata=zeros(tfmaxtimepts,maxfreqs,maxconds,maxpairs,ns);
phidata=zeros(tfmaxtimepts,maxfreqs,maxconds,maxpairs,ns);
cohstats=zeros(tfmaxtimepts,maxfreqs,maxconds,maxpairs,2);

%initialize a bad subject counter(?)
bad_s=zeros(ns,1);
cond_start = 1;

for o = 1:n_opt

n_conds = length(opt(o).case_vec);
cond_inds = cond_start:cond_start+n_conds-1;
cond_start = cond_start + n_conds;

%load in data
s_valid=0;
for s_attempt=1:ns
for datatype=1:length(datatypes)
load(matlist_all{s_attempt,o},datatypes{datatype});
end
%if exist('n_trials','var')
%    if ~( length(n_trials)==maxconds )
%        bad_s(s_attempt)=1;
%        continue
%    end
%end
s_valid=s_valid+1;
if exist('behav_data','var')
    behdata(s_valid)=behav_data;
end
if exist('coh','var')
if iscell(coh)
    cohdata(:,:,cond_inds,:,s_valid)=abs(coh{1}(ttfimes2import,:,:,:));
    cohstats(:,:,cond_inds,:,s_valid)=coh{2}(tftimes2import,:,:,:);
else
    cohdata(:,:,cond_inds,:,s_valid)=abs(coh(tftimes2import,:,:,:));
    phidata(:,:,cond_inds,:,s_valid)=angle(coh(tftimes2import,:,:,:));
end
end
if exist('n_trials','var')
n_trials_all(cond_inds,s_valid)=n_trials;
end
if exist('erp','var')
erpdata(:,:,cond_inds,s_valid)=erp(times2import,1:maxchans,:);
end
%if exist('wave_evk','var')
%wave_evkdata(:,:,:,:,s_valid)=wave_evk(tftimes2import,1:maxchans,:,:);
%end
if exist('wave_evknorm','var')
wave_evknormdata(:,:,:,cond_inds,s_valid)=abs(wave_evknorm(tftimes2import,1:maxchans,:,:));
end
%if exist('wave_tot','var')
%wave_totdata(:,:,:,cond_inds,s_valid)=wave_tot(tftimes2import,1:maxchans,:,:).^2; %square to give power (uV^2)
%end
if exist('wave_totpow','var')
wave_totpowdata(:,:,:,cond_inds,s_valid)=wave_totpow(tftimes2import,1:maxchans,:,:);
end
fprintf(num2str(s_attempt))
if mod(s_attempt,20)==0
    fprintf('\n')
end
end
%if nargout>=8
%    itcdata=abs(wave_evkdata)./wave_totdata;
%end
fprintf('\n')
end

%pack important import vars into struct
imp=v2struct(bad_s,maxchans,maxconds,maxfreqs,maxpairs,tfmaxtimepts,s_valid, ...
    tftimerate, erpmaxtimepts, erptimerate);

end