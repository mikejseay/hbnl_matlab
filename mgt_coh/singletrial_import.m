%% define paths for essential files

if true

singletrial_path='/active_projects/mike/pcheck_singletrial/ern_9_b1_a0000743_256_cnt_singletrial.mat';
trialmat_path='/active_projects/mike/fmri_phase4_ern_cleanCSD_newconds/ern_9_b1_a0000743_256_cnt.h1.61c_20f.mat';
h1_path='/processed_data/cnt-h1-files/ern/256/a-subjects/ns32-64/ern_9_b1_a0000743_256_cnt.h1';

% load single trial data and take absolute value

load(singletrial_path) %takes a while
st_tfdata=abs(st_tfdata);

end

%% load the trial indexing matrix and the full event table

load(trialmat_path,'trial_mat')
h1_struct = read_hdf1_behdata(h1_path);
etable=h1_getbehav(h1_struct,90);
for ev=1:size(etable,1)
    if ~isnan(etable.rt(ev))
        etable.rt(ev-1)=etable.rt(ev);
    end
end
etable=etable(etable.response_code==0,:);
etable=etable(2:end-1,:);

%%

resps=logical(sum(trial_mat(:,1:4),2));
RT=etable.rt(logical(resps));

%% heck, do some regressions on the identity of trials!!

n_conds=size(trial_mat,2);
[n_timepts,n_freqs,n_chans,n_trials]=size(st_tfdata);
n_pairs=size(st_cohdata,3);

beta=zeros(2,n_chans,n_timepts,n_freqs);
stats=zeros(4,n_chans,n_timepts,n_freqs);

for chan=1:n_chans
for timept=1:n_timepts
for freq=1:n_freqs
    [beta(:,chan,timept,freq), ~, ~, ~, stats(:,chan,timept,freq)]= ...
        regress(squeeze(st_tfdata(timept,freq,chan,resps)),...
        [ones(size(RT)) RT]);
end
end
end
% this ended up taking quite a long time

%% plot

for chan=[7 16 25]
    
    figure;
    imagesc(flipud(squeeze(beta(2,chan,:,:)))')
    
end