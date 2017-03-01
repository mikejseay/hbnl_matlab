%% condition means of TF results

% takes approximately 700 seconds per subject (~12 mins)

%load h1 file list
load('/export/home/mike/matlab/origin/mgt_coh/fmri_ern_h1_files.mat')
%params with all new steps
load('/export/home/mike/matlab/batch/ern_fmri_pchecknew.mat')
%

%%

for i=setdiff(15:length(h1_list),50)
[~,name,~]=fileparts(h1_list{i});
in_file=[param_struct.raw_path,'/',name,'_cleanraw.mat'];
load(in_file);

min_cycles=2.5;
cycle_linfactor=0.5;
maxfreq=30;
freq_range=[2 30];



[n_samps,n_epochs,n_chans]=size(trial_data);
[~,n_conds]=size(trial_mat);
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
data=reshape(trial_data(:,:,1),[1 n_samps n_epochs]);
[~, ~, ~, times, freqs, ~, ~, ~] = ...
            newtimef(data,n_samps,[ms_start,ms_end],fs,'cycles',[min_cycles cycle_linfactor], ...
            'maxfreq',maxfreq,'freqs',freq_range,'plotersp','off','plotitc','off');
n_times=length(times);
n_freqs=length(freqs);
end
%n_times=200;
%n_freqs=29;

erspdata=zeros(n_times,n_freqs,n_chans,n_conds);
itcdata=zeros(n_times,n_freqs,n_chans,n_conds);
cohdata=zeros(n_times,n_freqs,n_pairs,n_conds);

tic
for chan=1:n_chans
    for cond=1:n_conds
        data=reshape(trial_data(:,trial_mat(:,cond),chan),[1 n_samps sum(trial_mat(:,cond))]);
        %figure;
        %[ersp, itc, powbase, times, freqs, erspboot, itcboot, tfdata]
        [ersp, itc] = ...
            newtimef(data,n_samps,[ms_start,ms_end],fs,'cycles',[min_cycles cycle_linfactor], ...
            'maxfreq',maxfreq,'freqs',freq_range,'plotersp','off','plotitc','off');
        
        %assign the data into the matrix
        erspdata(:,:,chan,cond)=ersp';
        itcdata(:,:,chan,cond)=itc';
    end
end


for pair=1:n_pairs
    for cond=1:n_conds
        data1=reshape(trial_data(:,trial_mat(:,cond),pairs(pair,1)),[1 n_samps sum(trial_mat(:,cond))]);
        data2=reshape(trial_data(:,trial_mat(:,cond),pairs(pair,2)),[1 n_samps sum(trial_mat(:,cond))]);
        %[coh,mcoh,timesout,freqsout,cohboot,cohangles,allcoher,alltfX,alltfY]
        coh =...
            newcrossf(data1,data2,n_samps,[ms_start,ms_end],fs,[min_cycles cycle_linfactor], ...
            'maxfreq',maxfreq,'freqs',freq_range,'plotamp','off','plotphase','off');
        %assign the data into the matrix
        cohdata(:,:,pair,cond)=coh';
    end
end
toc

clear chan cond pair data data1 data2 ersp itc coh


tfparams=v2struct(min_cycles,cycle_linfactor,fs,freq_range,maxfreq, ...
    t_ms,n_chans,n_conds,n_pairs,n_samps,n_epochs,pairs, ...
    n_times,n_freqs,times,freqs);
clear min_cycles cycle_linfactor fs freq_range maxfreq ...
    t_ms n_chans n_conds n_pairs n_samps n_epochs pairs ...
    n_times n_freqs times freqs ms_end ms_start


Y=v2struct(behav_data,cohdata,erspdata,itcdata,mean_data,name,param_struct, ...
    tfparams);

out_file=[param_struct.condmean_path,'/',name,'_condmean.mat'];
save(out_file,'-struct','Y','-v6');
clear behav_data cohdata erspdata itcdata mean_data name tfparams
end
