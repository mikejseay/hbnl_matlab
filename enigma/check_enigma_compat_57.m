% check_enigma_compat_57.m - For 57 subjects, compares the enigma
% calculation results starting from two different original files: those
% that use MassComp calculations and those that have more relatively "raw"
% data from Neuroscan.

%load the 57 subjects info
load('enigma_57subs.mat')

%index the channels/names

%mc_chans=[20, 1, 16];
%ns_chans=[16, 31, 30];
%chan_names={'CZ','O1','O2'};

mc_chans=[18, 20, 21];
ns_chans=[7, 16, 25];
chan_names={'FZ','CZ','PZ'};

%choose an alpha range
a_range=[15:25];
t_range=[7:11];

%initialize a dummy var
good_s=0;

%max power val for display
maxpow_a=40;
maxpow_t=20;

%%

%preallocate variables
peak_a_mc=zeros(3,1); peak_a_mc_ind=zeros(3,1);
peak_a_ns=zeros(3,1); peak_a_ns_ind=zeros(3,1);
peak_t_mc=zeros(3,1); peak_t_mc_ind=zeros(3,1);
peak_t_ns=zeros(3,1); peak_t_ns_ind=zeros(3,1);
alphapeaks=zeros(2,55,3); thetapeaks=zeros(2,55,3);
PFMat_mc=zeros(55,3); PFMat_ns=zeros(55,3);
PowFMat_mc=zeros(55,5,3); PowFMat_ns=zeros(55,5,3);

%for each
for s=1:length(subjectid_table)

%load MassComp data, do enigma analysis, and save output as outputvar1
mc_out=enigma_analysis_mc(mc_files{s},mc_chans);

%load NeuroScan data, do enigma analysis, and save output as outputvar2
ns_out=enigma_analysis_ns(ns_files{s},ns_chans);

%if bad, skip
if ~isstruct(mc_out) || ~isstruct(ns_out)
    continue
end

%if we reached this step, it was a good subject so increment the counter
good_s=good_s+1;

%pull out peak alpha/theta values
for chan=1:size(mc_out.PW,2)
    %alpha
    [tpeak_a_mc,tpeak_a_mc_ind]=findpeaks(mc_out.PW(a_range,chan),'SORTSTR','descend');
    if isempty(tpeak_a_mc)
        [tpeak_a_mc,tpeak_a_mc_ind]=max(mc_out.PW(a_range,chan));
    end
    
    [tpeak_a_ns,tpeak_a_ns_ind]=findpeaks(ns_out.PW(a_range,chan),'SORTSTR','descend');
    if isempty(tpeak_a_ns)
        [tpeak_a_ns,tpeak_a_ns_ind]=max(ns_out.PW(a_range,chan));
    end
    peak_a_mc(chan)=tpeak_a_mc(1); peak_a_mc_ind(chan)=tpeak_a_mc_ind(1);
    peak_a_ns(chan)=tpeak_a_ns(1); peak_a_ns_ind(chan)=tpeak_a_ns_ind(1);
    %theta
    [tpeak_t_mc,tpeak_t_mc_ind]=findpeaks(mc_out.PW(t_range,chan),'SORTSTR','descend');
    if isempty(tpeak_t_mc)
        [tpeak_t_mc,tpeak_t_mc_ind]=max(mc_out.PW(t_range,chan));
    end
    
    [tpeak_t_ns,tpeak_t_ns_ind]=findpeaks(ns_out.PW(t_range,chan),'SORTSTR','descend');
    if isempty(tpeak_t_ns)
        [tpeak_t_ns,tpeak_t_ns_ind]=max(ns_out.PW(t_range,chan));
    end
    peak_t_mc(chan)=tpeak_t_mc(1); peak_t_mc_ind(chan)=tpeak_t_mc_ind(1);
    peak_t_ns(chan)=tpeak_t_ns(1); peak_t_ns_ind(chan)=tpeak_t_ns_ind(1);
end

%plot the spectrums side by side, highlighting peaks
%figure
%subplot(1,2,1); plot(mc_out.PW); grid on; axis([0 30 0 maxpow_a]);
%hold on; plot(a_range(1)+peak_a_mc_ind-1,peak_a_mc,'*'); hold off;
%hold on; plot(t_range(1)+peak_t_mc_ind-1,peak_t_mc,'m*'); hold off;
%subplot(1,2,2); plot(ns_out.PW); grid on; axis([0 30 0 maxpow_a]);
%hold on; plot(a_range(1)+peak_a_ns_ind-1,peak_a_ns,'*'); hold off;
%hold on; plot(t_range(1)+peak_t_ns_ind-1,peak_t_ns,'m*'); hold off;
%legend(chan_names)

%save the alpha/theta peak values for each
alphapeaks(1,good_s,:)=peak_a_mc;
alphapeaks(2,good_s,:)=peak_a_ns;
thetapeaks(1,good_s,:)=peak_t_mc;
thetapeaks(2,good_s,:)=peak_t_ns;

%save the PowF values for each
PowFMat_mc(good_s,:,:)=mc_out.PowF;
PowFMat_ns(good_s,:,:)=ns_out.PowF;

%save the PF values for each
PFMat_mc(good_s,:)=mc_out.PF;
PFMat_ns(good_s,:)=ns_out.PF;

end

%% check peak alpha/theta power differences across all subjects

%scatterplot the alpha peaks against each other for each channel, fit a
%line and save the slope/intercept
p=zeros(2,3,2);
figure(50)
for chan=1:3
    hold on; scatter(alphapeaks(1,:,chan),alphapeaks(2,:,chan));
    p(1,chan,:)=polyfit(alphapeaks(1,:,chan),alphapeaks(2,:,chan),1);
    hold on; plot(linspace(1,maxpow_a,maxpow_a),p(1,chan,1)*linspace(1,maxpow_a,maxpow_a)+p(1,chan,2))
end
axis([0 maxpow_a 0 maxpow_a])
hold on; plot(linspace(1,maxpow_a,maxpow_a),linspace(1,maxpow_a,maxpow_a),'k--')
hold off
%legend(chan_names)
title('Alpha')
xlabel('MassComp Alpha Power')
ylabel('NeuroScan Alpha Power')

%scatterplot the theta peaks against each other
figure(51)
for chan=1:3
    hold on; scatter(thetapeaks(1,:,chan),thetapeaks(2,:,chan));
    p(2,chan,:)=polyfit(thetapeaks(1,:,chan),thetapeaks(2,:,chan),1);
    hold on; plot(linspace(1,maxpow_t,maxpow_t),p(2,chan,1)*linspace(1,maxpow_t,maxpow_t)+p(2,chan,2))
end
axis([0 maxpow_t 0 maxpow_t])
hold on; plot(linspace(1,maxpow_t,maxpow_t),linspace(1,maxpow_t,maxpow_t),'k--')
hold off
%legend(chan_names)
title('Theta')
xlabel('MassComp Theta Power')
ylabel('NeuroScan Theta Power')

%% check power differences at each channel

%fit lines across channels for each subject
alpha_p_s=zeros(2,size(alphapeaks,2),2);
theta_p_s=zeros(2,size(thetapeaks,2),2);
for s=1:size(alphapeaks,2)
    alpha_p_s(1,s,:)=polyfit(1:3,squeeze(alphapeaks(1,s,:))',1); alpha_p_s(2,s,:)=polyfit(1:3,squeeze(alphapeaks(2,s,:))',1);
    theta_p_s(1,s,:)=polyfit(1:3,squeeze(thetapeaks(1,s,:))',1); theta_p_s(2,s,:)=polyfit(1:3,squeeze(thetapeaks(2,s,:))',1); 
end

%scatterplot the slope values
figure(52)
scatter(alpha_p_s(1,:,1),alpha_p_s(2,:,1))
xlabel('MassComp - Alpha')
ylabel('NeuroScan - Alpha')
title('Slope')
axis([0 25 0 25]); hold on;
plot(linspace(0,25,25),linspace(0,25,25),'k--'); grid on; hold off;

figure(53)
scatter(theta_p_s(1,:,1),theta_p_s(2,:,1))
xlabel('MassComp - Theta')
ylabel('NeuroScan - Theta')
title('Slope')
axis([0 3 0 3]); hold on;
plot(linspace(0,3,3),linspace(0,3,3),'k--'); grid on; hold off;

%scatterplot the intercept values
figure(54)
scatter(alpha_p_s(1,:,2),alpha_p_s(2,:,2))
xlabel('MassComp - Alpha')
ylabel('NeuroScan - Alpha')
title('Intercept')
axis([0 25 0 25]); hold on;
plot(linspace(0,25,25),linspace(0,25,25),'k--'); grid on; hold off;

figure(55)
scatter(theta_p_s(1,:,2),theta_p_s(2,:,2))
xlabel('MassComp - Theta')
ylabel('NeuroScan - Theta')
title('Intercept')
axis([0 25 0 25]); hold on;
plot(linspace(0,25,25),linspace(0,25,25),'k--'); grid on; hold off;


%% check differences in the measured PowF variable by the enigma script

freqlabels={'Delta','Theta','Alpha','Beta','Full Spectrum'};

figure(56)
xlabel('MassComp')
ylabel('NeuroScan')
for freq=1:5
    for chan=1:3
        subplot(1,5,freq)
        scatter(PowFMat_mc(:,freq,chan),PowFMat_ns(:,freq,chan))
        hold on
        p_powf(chan,freq,:)=polyfit(PowFMat_mc(:,freq,chan),PowFMat_ns(:,freq,chan),1);
        plot(linspace(-4,4,9),linspace(-4,4,9)*p_powf(chan,freq,1)+p_powf(chan,freq,2))
        hold on
    end
    plot(linspace(-4,4,9),linspace(-4,4,9),'k--')
    hold off
    axis([-4 4 -4 4])
    title(freqlabels{freq})
    grid on
end
hold off


%% check differences in the measured PF variable by enigma

figure(57)
for chan=1:3
    scatter(PFMat_mc(:,chan),PFMat_ns(:,chan))
    hold on
    p_pf(chan,:)=polyfit(PFMat_mc(:,chan),PFMat_ns(:,chan),1);
    plot(linspace(8.8,10.8,3),linspace(8.8,10.8,3)*p_pf(chan,1)+p_pf(chan,2))
    hold on
end
plot(linspace(8.8,10.8,3),linspace(8.8,10.8,3),'k--')
hold off