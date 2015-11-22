function [mat_list, imp, n_trials_all, erpdata, wave_evkdata, wave_totdata, cohdata, itcdata] = ...
    coh_import_stage(opt, demogsfile, time_downsamp)
% Imports computed EEG measures into the MATLAB workspace for freestyle
% plotting and exploration.

if nargin < 3
    time_downsamp=2; % downsampling by half is default
end


% make a list of each subject's .mat file
if nargin < 2 % if no demogsfile was specified, use this
    mat_list=coh_updatemats(opt.outpath);
else
    mat_list=coh_updatemats(opt.outpath,demogsfile);
end


ns=length(mat_list);

%size of data dimensions, taking into account a time-downsampling factor
maxtimepts=round((opt.n_samps)/time_downsamp);
times2import=1:time_downsamp:maxtimepts*time_downsamp;
timerate=opt.rate/time_downsamp;

maxfreqs=length(opt.wavelet_scales);
maxconds=length(opt.case_vec);
maxchans=length(opt.chan_vec);
maxpairs=size(opt.coherence_pairs,1);

%indicate mats to load
n_datatypes=min(nargout-2, 5);
datatype_names={'n_trials','erp','wave_evk','wave_tot','coh'};
datatypes={datatype_names{1:n_datatypes}};

%pre-allocate vars
n_trials_all=zeros(maxconds,ns);
erpdata=zeros(maxtimepts,maxchans,maxconds,ns);
wave_evkdata=zeros(maxtimepts,maxchans,maxfreqs,maxconds,ns);
wave_totdata=zeros(maxtimepts,maxchans,maxfreqs,maxconds,ns);
cohdata=zeros(maxtimepts,maxfreqs,maxconds,maxpairs,ns);
cohstats=zeros(maxtimepts,maxfreqs,maxconds,maxpairs,2);

%initialize a bad subject counter(?)
bad_s=zeros(ns,1);

%load in data
s_valid=0;
for s_attempt=1:ns
for datatype=1:length(datatypes)
load(mat_list{s_attempt},datatypes{datatype});
end
if exist('n_trials','var')
    if ~( length(n_trials)==maxconds )
        bad_s(s_attempt)=1;
        continue
    end
end
s_valid=s_valid+1;
if exist('behav_data','var')
    behdata(s_valid)=behav_data;
end
if exist('coh','var')
if iscell(coh)
    cohdata(:,:,:,:,s_valid)=abs(coh{1}(times2import,:,:,:));
    cohstats(:,:,:,:,s_valid)=coh{2}(times2import,:,:,:);
else
    cohdata(:,:,:,:,s_valid)=abs(coh(times2import,:,:,:));
end
end
if exist('n_trials','var')
n_trials_all(:,s_valid)=n_trials;
end
if exist('erp','var')
erpdata(:,:,:,s_valid)=erp(times2import,1:maxchans,:);
end
if exist('wave_evk','var')
wave_evkdata(:,:,:,:,s_valid)=wave_evk(times2import,1:maxchans,:,:);
wave_totdata(:,:,:,:,s_valid)=wave_tot(times2import,1:maxchans,:,:);
end
fprintf(num2str(s_attempt))
if mod(s_attempt,20)==0
    fprintf('\n')
end
end
if nargout==8
    itcdata=abs(wave_evkdata)./wave_totdata;
end
fprintf('\n')

%pack important import vars into struct
imp=v2struct(bad_s,maxchans,maxconds,maxfreqs,maxpairs,maxtimepts,s_valid, ...
    timerate);

%load('/export/home/mike/matlab/origin/fmri/fmri_ern_beh.mat')
if exist('behdata','var')
    beh_import
end

end