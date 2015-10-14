% set certain plotting parameters for coh_analysis

cond_diff={[2],[1]}; %[6],[7]};
scl.cond_label{imp.maxconds+1}= ...
    [scl.cond_label{cond_diff{1}},' - ',scl.cond_label{cond_diff{2}}];
%scl.cond_label{imp.maxconds+1}='cond diff';



%f_start_hz=[2.9,3.6,6.4,8]; %lower edges
%f_end_hz=[3.2,5.3,8,12.8]; %upper edges
%f_indiv_hz=[3.2, 4.6, 7.1, 11,  18];        %indiv bands

chosen_s=1;
coh_absmin=0;
coh_absmax=1;

%[9 7 8 17 16 18 23 25 24];
% for eros [7 8 16 18 25];
chosen_chan=[9 7 8 17 16 18 23 25 24];
chosen_cond=1;


if ~exist('tfparams','var')
f_start_hz=[2.9, 3.6, 4.6, 6.4, 8,    4];     %lower edges
f_end_hz=[  3.2, 5.3, 6.4, 8,   12.8, 5.3];  %upper edges
f_indiv_hz=[3.2, 4.6, 5.3, 7.1, 10.7, 4.6]; %indiv bands
else
f_start_hz=[2, 3, 5, 8,  13];   %lower edges
f_end_hz=[  2, 5, 7, 12, 24];   %upper edges
f_indiv_hz=[2, 4, 6, 10, 18];   %indiv bands
end
chosen_freq=[6 10 13 16 19];

chosen_g=[10 11];
chosen_p=[14 12 4 9];

%t_start_ms=[0 100 200 300 400 600];
%t_end_ms=[100 200 300 400 600 800];
t_win_start=0;
t_win_end=500;
t_win_space=100;
t_start_ms=[t_win_start:t_win_space:t_win_end];
t_end_ms=t_start_ms+t_win_space;

plotn_chan=1:9; %[2 5 8]; %[1 2 3 4 6 8]; %3
plotn_cond=1:imp.maxconds+1;
plotn_f=[1 3 5];
plotn_g=1:2;
plotn_p=[1:4];
plotn_t=1:length(t_start_ms);


maxwin=length(t_start_ms);
%hist_ax=[-1 1 0 round(size(s_inds_g,1)/6)];
hist_ax=[-.8 .8 0 10];
%hist_nbin=round(size(s_inds_g,1)/6);
hist_nbin=10;
p_loc=[hist_ax(1)/2 hist_ax(4)*3/4];
cmap_bone=flipud(bone(256));
pmkmp_scheme='cubicl'; %'linlhot' is also good
cmap=pmkmp(256,pmkmp_scheme);
%cmap(1:48,:)=[];
pmkmp_scheme_line='linlhot'; %'linlhot' is also good
cmap_line=pmkmp(256,pmkmp_scheme_line);
figdum_init=1;
figdum=figdum_init;
n_contour=8;

chosen_topochan=1:61;

chosen_coh_limit=[0.25 0.4];
chosen_coh_limit_diff=[-.1 .1];

%subplot structure (based on conditions)
sp_d=numSubplots(length(plotn_cond));
%sp_d=[1 length(plotn_cond)];

%pack into a struct
pp=v2struct(cond_diff, chosen_g, chosen_chan, f_start_hz, f_end_hz, chosen_freq, ...
    chosen_cond, chosen_s, coh_absmin, coh_absmax, t_start_ms, t_end_ms, maxwin, hist_ax, ...
    hist_nbin, p_loc, cmap_bone, cmap, cmap_line, pmkmp_scheme_diff, ...
    pmkmp_scheme_line, figdum_init, figdum, n_contour, ...
    plotn_chan, plotn_cond, plotn_g, plotn_f, plotn_t, chosen_topochan, ...
    chosen_coh_limit, chosen_coh_limit_diff, f_indiv_hz, chosen_p, plotn_p);
clear cond_diff chosen_g chosen_chan f_start_hz f_end_hz chosen_freq ...
    chosen_cond chosen_s coh_absmin coh_absmax t_start_ms t_end_ms maxwin hist_ax ...
    hist_nbin p_loc pmkmp_scheme cmap figdum figdum_init n_contour ...
    plotn_chan plotn_cond plotn_g plotn_f plotn_t chosen_topochan ...
    chosen_coh_limit chosen_coh_limit_diff f_indiv_hz chosen_p plotn_p ...
    cmap_line pmkmp_scheme_line

clear t_win_start t_win_end t_win_space