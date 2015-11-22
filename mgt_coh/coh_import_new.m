%declare number of subjects to look at
ns=length(mat_list);

%figure out the data dimensions
load(mat_list{1})
[maxtimepts,maxfreqs,maxpairs,maxconds]=size(cohdata);
maxchans=size(erspdata,3);

%indicate variables to load from .mat
datatypes={'cohdata','erspdata','itcdata'};

%pre-allocate vars
cohdata=zeros(maxtimepts,maxfreqs,maxconds,maxpairs,2);
%n_trials_all=zeros(maxconds,2);
wave_totdata=zeros(maxtimepts,maxchans,maxfreqs,maxconds,2);
wave_evkdata=zeros(maxtimepts,maxchans,maxfreqs,maxconds,2);

%initialize a bad subject counter(?)
bad_s=zeros(ns,1);

%add prefix if being run on pc
pc_netdrive_prefix='';
if ispc
    for s_attempt=1:ns
        mat_list{s_attempt}=[pc_netdrive_prefix,mat_list{s_attempt}];
    end
end

%load in data
s_valid=0;
for s_attempt=1:ns
cohdata_load=load(mat_list{s_attempt},'cohdata');
cohdata_load=cohdata_load.cohdata;
erspdata_load=load(mat_list{s_attempt},'erspdata');
erspdata_load=erspdata_load.erspdata;
itcdata_load=load(mat_list{s_attempt},'itcdata');
itcdata_load=itcdata_load.itcdata;
load(mat_list{s_attempt},'behav_data');
%if ~( length(n_trials)==maxconds )
%    bad_s(s_attempt)=1;
%    continue
%end
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
cohdata(:,:,:,:,s_valid)=permute(cohdata_load,[1 2 4 3]);
wave_totdata(:,:,:,:,s_valid)=permute(erspdata_load,[1 3 2 4]);
wave_evkdata(:,:,:,:,s_valid)=permute(itcdata_load,[1 3 2 4]);
end
%flip the frequency dimension for plots
%cohdata=flip(cohdata,2);
%wave_totdata=flip(wave_totdata,2);
%itcdata=flip(itcdata,2);
itcdata=abs(wave_evkdata);

%pack important import vars into struct
imp=v2struct(bad_s,maxchans,maxconds,maxfreqs,maxpairs,maxtimepts,s_valid);
clear bad_s maxchans maxconds maxfreqs maxpairs maxtimepts s_valid

clear ns coh_stats_present coh_results mean_data n_trials wavelet_evk ...
    wavelet_tot s_attempt datatype datatypes pc_netdrive_prefix ...
    behav_data ans cohdata_load erspdata_load itcdata_load

if exist('behdata','var')
    beh_import
end