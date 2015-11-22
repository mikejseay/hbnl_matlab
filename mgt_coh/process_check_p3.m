%% single trial TF results

% takes approximately 1600 seconds (27 mins) / subject

%load h1 file list
load('/export/home/mike/matlab/origin/mgt_coh/fmri_ern_h1_files.mat')
%params with all new steps
load('/export/home/mike/matlab/batch/ern_fmri_pchecknew.mat')
%
[~,name,~]=fileparts(h1_list{1});
in_file=[param_struct.raw_path,'/',name,'_cleanraw.mat'];
load(in_file);

min_cycles=2.5;
cycle_linfactor=0.5;
maxfreq=30;
freq_range=[2 30];

%%

[n_samps,n_trials,n_chans]=size(trial_data);
pairs=param_struct.coherence_pairs;
n_pairs=length(pairs);

fs=param_struct.rate;
ms_start=-param_struct.prestim_ms;
ms_end=ms_start + ( n_samps * (1000 / fs) );
t_ms=linspace(ms_start, ms_end, n_samps+1);
t_ms=t_ms(1:n_samps);

%do a quick newtimef with the desired parameters in order to determine the
%matrix size for the outputs

if true
data=reshape(trial_data(:,:,1),[1 n_samps n_trials]);
[~, ~, ~, times, freqs] = ...
            newtimef(data,n_samps,[ms_start,ms_end],fs,'cycles',[min_cycles cycle_linfactor], ...
            'maxfreq',maxfreq,'freqs',freq_range,'plotersp','off','plotitc','off');
n_times=length(times);
n_freqs=length(freqs);
end
%n_times=200;
%n_freqs=29;

st_tfdata=zeros(n_times,n_freqs,n_chans,n_trials);
st_cohdata=zeros(n_times,n_freqs,n_pairs,n_trials);

tic
for chan=1:n_chans
    data=reshape(trial_data(:,:,chan),[1 n_samps n_trials]);
    %figure;
    %[ersp, itc, powbase, times, freqs, erspboot, itcboot, st_tfdata]
    [~, ~, ~, ~, ~, ~, ~, chan_tfdata] = ...
        newtimef(data,n_samps,[ms_start,ms_end],fs,'cycles',[min_cycles cycle_linfactor], ...
        'maxfreq',maxfreq,'freqs',freq_range,'plotersp','off','plotitc','off');
    chan_tfdata=permute(chan_tfdata,[2 1 3]);
    %assign the data into the matrix
    st_tfdata(:,:,chan,:)=chan_tfdata;
end

for pair=1:n_pairs
    data1=reshape(trial_data(:,:,pairs(pair,1)),[1 n_samps n_trials]);
    data2=reshape(trial_data(:,:,pairs(pair,2)),[1 n_samps n_trials]);
    %[coh,mcoh,timesout,freqsout,cohboot,cohangles,allcoher,alltfX,alltfY]
    [~,~,~,~,~,~,allcoher] =...
        newcrossf(data1,data2,n_samps,[ms_start,ms_end],fs,[min_cycles cycle_linfactor], ...
        'maxfreq',maxfreq,'freqs',freq_range,'plotamp','off','plotphase','off');
    %assign the data into the matrix
    allcoher=permute(allcoher,[2 1 3]);
    st_cohdata(:,:,pair,:)=allcoher;
end
toc

clear chan pair data data1 data2 chan_tfdata allcoher

%%

tfparams=v2struct(min_cycles,cycle_linfactor,fs,freq_range,maxfreq, ...
    t_ms,n_chans,n_pairs,n_samps,n_trials,pairs, ...
    n_times,n_freqs,times,freqs);
clear min_cycles cycle_linfactor fs freq_range maxfreq ...
    t_ms n_chans n_conds n_pairs n_samps n_trials pairs ...
    n_times n_freqs times freqs ms_end ms_start

%%

Y=v2struct(behav_data,st_tfdata,st_cohdata,name,param_struct, ...
    tfparams);

out_file=[param_struct.singletrial_path,'/',name,'_singletrial.mat'];
save(out_file,'-struct','Y','-v6');
clear behav_data st_tfdata st_cohdata name param_struct ...
    tfparams