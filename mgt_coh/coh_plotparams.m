function [scl, pp] = coh_plotparams(imp, scl, cond_diff)
% set certain plotting parameters for coh_analysis

%cond_diff={[1 2],[3 4]}; %[6],[7]};
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
%chosen_chan=[7 8 9 16 17 18 23 24 25 56 24 55];
%chosen_chan=[58 27 26 25 56 7 42 41 15 14 16 23 24 38];
chosen_chan=[9 7 8 23 25 24 31 58 30];
chosen_cond=1;


if ~exist('tfparams','var')
%f_start_hz=[2.9, 3.6, 4.6, 6.4, 8,    4];     %lower edges
%f_end_hz=[  3.2, 5.3, 6.4, 8,   12.8, 5.3];  %upper edges
%f_indiv_hz=[3.2, 4.6, 5.3, 7.1, 10.7, 4.6]; %indiv bands

%f_start_hz=[2,   3.8, 5.5, 8,    14,   27  ];     %lower edges
%f_end_hz=[  3.5, 5,   7.3, 13,   24,   32  ];  %upper edges
%f_indiv_hz=[2.6, 4.2, 6.1, 9.7,  18.3, 28.4 ]; %indiv bands

%f_start_hz=[2,   3.5, 4.2, 5.5,  8,    14,   27, 3, 5, 4];     %lower edges
%f_end_hz=[  3.2, 4.2, 5,   7.3,  13,   24,   32, 5, 7, 6];  %upper edges
%f_indiv_hz=[2.6, 3.9, 4.6, 6.1,  10,  18, 28,    4, 6, 5]; %indiv bands

f_start_hz=[2,   4,   8,    14, 26, 3, 5];     %lower edges
f_end_hz=[  3,   7,   13,   25, 32, 5, 7];  %upper edges
f_indiv_hz=[2.5, 5.5, 10.5, 19, 29, 4, 6]; %indiv bands

else
f_start_hz=[2, 3, 5, 9,  13];   %lower edges
f_end_hz=[  2, 5, 7, 12, 24];   %upper edges
f_indiv_hz=[2, 4, 6, 10, 18];   %indiv bands
end
chosen_freq=[6 10 13 16 19];

chosen_g=[2 3]; %[1]; %
chosen_p=[14 12 4 9];

%t_start_ms=[0 100 200 300 400 600];
%t_end_ms=[100 200 300 400 600 800];

t_win_start=0;
t_win_end=500;
t_win_space=100;
t_start_ms=[t_win_start:t_win_space:t_win_end];
t_end_ms=t_start_ms+t_win_space;

%t_start_ms=[90 160 250 290 340 400 500]; %-600
%t_end_ms=[110 180 270 310 380 500 600]; %-200

%t_start_ms=[0 200 400 600];
%t_end_ms=[200 400 600 800];

t_start_b_ms=-500;
t_end_b_ms=-200;


plotn_chan=1:length(chosen_chan); %[1 3 5 6 7 12 13 14]; %[1 3 5]; %[2 5 8]; %[1 2 3 4 6 8]; %3
plotn_cond=1:imp.maxconds+1;
plotn_f=1:length(f_start_hz);
plotn_g=1:length(chosen_g);
plotn_p=[1:4];
plotn_t=1:length(t_start_ms);


maxwin=length(t_start_ms);
%hist_ax=[-1 1 0 round(size(s_inds_g,1)/6)];
hist_ax=[-.6 .6 0 100];
%hist_nbin=round(size(s_inds_g,1)/6);
hist_nbin=30;
p_loc=[hist_ax(1)/2 hist_ax(4)*3/4];
%cmap_bone=flipud(bone(256));
pmkmp_scheme='cubicl'; %'linlhot' is also good
cmap=pmkmp(256,pmkmp_scheme);
%cmap(1:48,:)=[];
pmkmp_scheme_line='linlhot'; %'linlhot' is also good
cmap_line=pmkmp(256,pmkmp_scheme_line);
figdum_init=1;
figdum=figdum_init;
n_contour=8;

chosen_topochan=1:61;
topo_elecs='off';

chosen_coh_limit=[0.25 0.4];
chosen_coh_limit_diff=[-.1 .1];

%subplot structure (based on conditions)
sp_d=numSubplots(length(plotn_cond));
%sp_d=[1 length(plotn_cond)];

%pack into a struct
pp=v2struct(cond_diff, chosen_g, chosen_chan, f_start_hz, f_end_hz, chosen_freq, ...
    chosen_cond, chosen_s, coh_absmin, coh_absmax, t_start_ms, t_end_ms, ...
    t_start_b_ms, t_end_b_ms, maxwin, hist_ax, ...
    hist_nbin, p_loc, cmap, cmap_line, ...
    pmkmp_scheme_line, figdum_init, figdum, n_contour, ...
    plotn_chan, plotn_cond, plotn_g, plotn_f, plotn_t, chosen_topochan, topo_elecs, ...
    chosen_coh_limit, chosen_coh_limit_diff, f_indiv_hz, chosen_p, plotn_p, ...
    sp_d);

end