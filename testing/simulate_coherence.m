%understand wavelets and phase calculation!

%determine basic aspects of data
n_channels=10;
n_timepts_per_trial=256;
n_trials=100;

%build sine wave time signal per trial
srate=256;
dt=1/srate;
stoptime=1;
t=(0:dt:stoptime-dt)';

%assign amplitudes, frequencies, and phase lags for 3 sine components
a1=1; a2=0.5; a3=2;
f1=4; f2=8; f3=16;
p1=0; p2=pi/2; p3=pi;

%create fake data as a sum of sinusoids of various frequencies
n_total_timepts=n_timepts_per_trial*n_trials;
data=zeros(n_total_timepts,n_channels);
for chan=1:n_channels
    for trial=1:n_trials
        trial_pts=(1+(trial-1)*n_timepts_per_trial):(trial*n_timepts_per_trial);
        s1=a1*cos(2*pi*f1*t+p1);
        s2=a2*cos(2*pi*f2*randn*t+p2);
        s3=a3*cos(2*pi*f3*t+p3*randn);
        data(trial_pts,chan)=s1+s2+s3;
    end
end

%create a scale of frequencies which are 2*srate/freq
scale=[14,16,18,20,24,28,32,36,40,48,56,64,72,80,96,112,128,144,160,176];
n_scales=length(scale);

%do the wavelet calculation on the continuous data
dataW = wavelet_calc(data,scale);

%reshape the data into trials
dataWR=reshape(dataW, [n_timepts_per_trial,n_trials,n_channels,n_scales]);

%calculate ITC
itc=abs(mean(dataWR,2))./mean(abs(dataWR),2);

%define pairs of electrodes to examine, the trial structure(none), 
pairs=[1 2; 3 4; 5 6; 7 8; 9 10];
n_pairs=size(pairs,1);
trial_mat=ones(100,1);
n_perms=100;
seed=10;

%initialize the matrices for the coherence calculations
coh=zeros(n_timepts_per_trial,n_scales,n_pairs);
coh_stats=zeros(n_timepts_per_trial,n_scales,n_pairs);
coh_b=zeros(n_timepts_per_trial,n_scales,n_pairs);

%for each pair calculate the pairwise coherence
for pair=1:n_pairs
    [coh(:,:,pair),coh_stats(:,:,pair),B]=pairwise_coherence_calc(dataWR(:,:,pairs(pair,:),:),trial_mat,n_perms,seed);
end
coh_results=abs(coh);