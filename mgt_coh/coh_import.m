%declare number of subjects to look at
ns=length(mat_list);

%down-sample time-points by half?
time_downsamp=true;

%size of data dimensions
if time_downsamp
    maxtimepts=round((param_struct.n_samps-8)/2);
    times2import=1:2:maxtimepts*2;
    param_struct.timerate=param_struct.rate/2;
else
    maxtimepts=param_struct.n_samps-8;
    times2import=1:maxtimepts;
    param_struct.timerate=param_struct.rate;
end
maxfreqs=length(param_struct.scale);
maxconds=length(param_struct.case_vec);
maxchans=length(param_struct.chan_vec);
maxpairs=size(param_struct.coherence_pairs,1);

%indicate mats to load
datatypes={'mean_data','n_trials','coh_results'}; %'wavelet_evk','wavelet_tot'}; %,'coh_results'}; %,'behav_data'};
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
if exist('coh_results','var')
if iscell(coh_results)
    cohdata(:,:,:,:,s_valid)=abs(coh_results{1}(times2import,:,:,:));
    cohstats(:,:,:,:,s_valid)=coh_results{2}(times2import,:,:,:);
else
    cohdata(:,:,:,:,s_valid)=abs(coh_results(times2import,:,:,:));
end
end
if exist('n_trials','var')
n_trials_all(:,s_valid)=n_trials;
end
erpdata(:,:,:,s_valid)=mean_data(times2import,1:maxchans,:);
if exist('wavelet_evk','var')
wave_evkdata(:,:,:,:,s_valid)=wavelet_evk(times2import,1:maxchans,:,:);
wave_totdata(:,:,:,:,s_valid)=wavelet_tot(times2import,1:maxchans,:,:);
end
fprintf(num2str(s_attempt))
if mod(s_attempt,20)==0
    fprintf('\n')
end
end
%itcdata=abs(wave_evkdata)./wave_totdata;
fprintf('\n')
%filter ERP data
%erpdata=filter_erp_datastruct(erpdata,256,16);

%pack important import vars into struct
imp=v2struct(bad_s,maxchans,maxconds,maxfreqs,maxpairs,maxtimepts,s_valid);
clear bad_s maxchans maxconds maxfreqs maxpairs maxtimepts s_valid

clear ns coh_stats_present coh_results mean_data n_trials wavelet_evk ...
    wavelet_tot s_attempt datatype datatypes pc_netdrive_prefix ...
    behav_data trial_mat etable

%load('/export/home/mike/matlab/origin/fmri/fmri_ern_beh.mat')
if exist('behdata','var')
    beh_import
end