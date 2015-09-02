% set certain plotting parameters for coh_analysis

cond_diff={[9],[5]};
scl.cond_label{imp.maxconds+1}= ...
    [scl.cond_label{cond_diff{1}},' - ',scl.cond_label{cond_diff{2}}];
%scl.cond_label{imp.maxconds+1}='cond diff';



chosen_g=[10 11];
chosen_chan=[7 16 25 21 22];
chosen_p=[4 11 14];

%f_start_hz=[2.9,3.6,6.4,8]; %lower edges
%f_end_hz=[3.2,5.3,8,12.8]; %upper edges
f_start_hz=[2.9, 3.6, 4.6, 6.4, 8];    %lower edges
f_end_hz=[  3.2, 5.3, 6.4, 8,   12.8]; %upper edges
f_indiv_hz=[3.2, 4.6, 5.3, 7.1,  10.7];        %indiv bands
%f_indiv_hz=[3.2, 4.6, 7.1, 11,  18];        %indiv bands

chosen_freq=[6 10 13 16 19];
chosen_cond=1;
chosen_s=1;
coh_absmin=0;
coh_absmax=1;

%t_start_ms=[0 100 200 300 400 600];
%t_end_ms=[100 200 300 400 600 800];
t_win_start=0;
t_win_end=400;
t_win_space=100;
t_start_ms=[t_win_start:t_win_space:t_win_end];
t_end_ms=t_start_ms+t_win_space;

maxwin=length(t_start_ms);
%hist_ax=[-1 1 0 round(size(s_inds_g,1)/6)];
hist_ax=[-.8 .8 0 10];
%hist_nbin=round(size(s_inds_g,1)/6);
hist_nbin=10;
p_loc=[hist_ax(1)/2 hist_ax(4)*3/4];
pmkmp_scheme='cubicl'; %'linlhot' is also good
cmap=pmkmp(256,pmkmp_scheme);
figdum_init=1;
figdum=figdum_init;
n_contour=9;

plotn_chan=2; %3
plotn_p=1:3;
plotn_cond=5:10;
plotn_g=1:2;
plotn_f=[1 3];
plotn_t=1:length(t_start_ms);

%subplot structure (based on conditions)
sp_d=numSubplots(length(plotn_cond));

chosen_topochan=1:61;

chosen_coh_limit=[0 1];
chosen_coh_limit_diff=[-.5 .5];

%pack into a struct
pp=v2struct(cond_diff, chosen_g, chosen_chan, f_start_hz, f_end_hz, chosen_freq, ...
    chosen_cond, chosen_s, coh_absmin, coh_absmax, t_start_ms, t_end_ms, maxwin, hist_ax, ...
    hist_nbin, p_loc, pmkmp_scheme, cmap, figdum_init, figdum, n_contour, ...
    plotn_chan, plotn_cond, plotn_g, plotn_f, plotn_t, chosen_topochan, ...
    chosen_coh_limit, chosen_coh_limit_diff, f_indiv_hz, chosen_p, plotn_p);
clear cond_diff chosen_g chosen_chan f_start_hz f_end_hz chosen_freq ...
    chosen_cond chosen_s coh_absmin coh_absmax t_start_ms t_end_ms maxwin hist_ax ...
    hist_nbin p_loc pmkmp_scheme cmap figdum figdum_init n_contour ...
    plotn_chan plotn_cond plotn_g plotn_f plotn_t chosen_topochan ...
    chosen_coh_limit chosen_coh_limit_diff f_indiv_hz chosen_p plotn_p

clear t_win_start t_win_end t_win_space