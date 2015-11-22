%declare number of subjects to look at
ns=length(mat_list);

%down-sample time-points by half?
time_downsamp=true;

%size of data dimensions
if time_downsamp
    maxtimepts=round((opt.n_samps)/2);
    times2import=1:2:maxtimepts*2;
    opt.timerate=opt.rate/2;
else
    maxtimepts=opt.n_samps;
    times2import=1:maxtimepts;
    opt.timerate=opt.rate;
end
maxfreqs=length(opt.wavelet_scales);
maxconds=length(opt.case_vec);
maxchans=length(opt.chan_vec);
maxpairs=size(opt.coherence_pairs,1);

%indicate mats to load
%datatypes={'n_trials','erp','wave_evk','wave_tot','coh'}; %'wavelet_evk','wavelet_tot'}; %,'coh_results'}; %,'behav_data'};
datatypes={'n_trials','mean_data','wave_evk','wave_tot','coh'};
%datatypes='''mean_data'',''n_trials'',''coh_results'',''wavelet_evk'',''wavelet_tot''';
%datatypes={'behav_data'};

%pre-allocate vars
cohdata=zeros(maxtimepts,maxfreqs,maxconds,maxpairs,ns);
%cohstats=zeros(maxtimepts,maxfreqs,maxconds,maxpairs,2);
erpdata=zeros(maxtimepts,maxchans,maxconds,ns);
n_trials_all=zeros(maxconds,ns);
%wave_evkdata=zeros(maxtimepts,maxchans,maxfreqs,maxconds,ns);
%wave_totdata=zeros(maxtimepts,maxchans,maxfreqs,maxconds,ns);

%initialize a bad subject counter(?)
bad_s=zeros(ns,1);

%add prefix if being run on pc
%pc_netdrive_prefix='B:';
%if ispc
%    for s_attempt=1:ns
%        mat_list{s_attempt}=[pc_netdrive_prefix,mat_list{s_attempt}];
%    end
%end

%load in data
s_valid=0;
for s_attempt=1:ns
for datatype=1:length(datatypes)
load(mat_list{s_attempt},datatypes{datatype});
end
%eval(['load(mat_list{s_attempt},',strjoin(datatypes,','),');']);
%eval(['load(mat_list{s_attempt},',datatypes,');']);
if exist('n_trials','var')
    if ~( length(n_trials)==maxconds )
        bad_s(s_attempt)=1;
        continue
    end
end
if exist('behav_data','var')
    if any(any(isnan(behav_data.respRT)))
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
itcdata=abs(wave_evkdata)./wave_totdata;
fprintf('\n')
%filter ERP data
%erpdata=filter_erp_datastruct(erpdata,256,16);

%pack important import vars into struct
imp=v2struct(bad_s,maxchans,maxconds,maxfreqs,maxpairs,maxtimepts,s_valid);

clear bad_s maxchans maxconds maxfreqs maxpairs maxtimepts s_valid ...
    ns coh_stats_present coh_results mean_data n_trials wavelet_evk ...
    wavelet_tot s_attempt datatype datatypes pc_netdrive_prefix ...
    behav_data trial_mat etable time_downsamp erp wave_evk wave_tot coh ...
    times2import

%load('/export/home/mike/matlab/origin/fmri/fmri_ern_beh.mat')
if exist('behdata','var')
    beh_import
end