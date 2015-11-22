function pp = coh_plotparams_gui (imp, scl, pp)
% set certain plotting parameters for coh_analysis with a GUI

cond_diff={[2],[1]};
scl.cond_label{imp.maxconds+1}= ...
    [scl.cond_label{cond_diff{1}},' - ',scl.cond_label{cond_diff{2}}];

coh_absmin=0;
coh_absmax=1;

chosen_s=1;
chosen_chan=[7 8 9 16 17 18 23 24 25];
chosen_cond=1;
chosen_freq=[6 10 13 16 19];
chosen_g=[10 11];
chosen_p=[14 12 4 9];
chosen_topochan=1:61;
chosen_coh_limit=[0.25 0.4];
chosen_coh_limit_diff=[-.1 .1];

           %1   %2   %3   %4    %5    %6    %7  %8 %9
f_start_hz=[2,   3.5, 4.2, 5.5,  8,    14,   27, 3, 5]; %lower edges
f_end_hz=[  3.2, 4.2, 5,   7.3,  13,   24,   32, 5, 7]; %upper edges
f_indiv_hz=[2.6, 3.9, 6.1, 4.6,  10,  18, 28,    4, 6]; %indiv bands

t_win_start=-100;
t_win_end=500;
t_win_space=100;
t_start_ms=[t_win_start:t_win_space:t_win_end];
t_end_ms=t_start_ms+t_win_space;
%t_start_ms=[90 160 250 290 340 400 500];
%t_end_ms=[110 180 270 310 380 500 600];
%t_start_ms=[90 200 500];
%t_end_ms=[180 400 600];

plotn_chan=[1 4 9];
plotn_cond=[1 2 3];
plotn_f=1:7;
plotn_g=1:length(chosen_g);
plotn_p=[1:4];
plotn_t=1:length(t_start_ms);

maxwin=length(t_start_ms);
hist_ax=[-.6 .6 0 100];
hist_nbin=30;
p_loc=[hist_ax(1)/2 hist_ax(4)*3/4];
pmkmp_scheme='cubicl'; %'linlhot' is also good
cmap=pmkmp(256,pmkmp_scheme);
pmkmp_scheme_line='linlhot';
cmap_line=pmkmp(256,pmkmp_scheme_line);
figdum_init=1;
figdum=figdum_init;
n_contour=8;

%subplot structure (based on conditions)
sp_d=numSubplots(length(plotn_cond));

%pack into a struct
pp=v2struct(cond_diff, chosen_g, chosen_chan, f_start_hz, f_end_hz, chosen_freq, ...
    chosen_cond, chosen_s, coh_absmin, coh_absmax, t_start_ms, t_end_ms, maxwin, hist_ax, ...
    hist_nbin, p_loc, cmap, cmap_line, ...
    pmkmp_scheme_line, figdum_init, figdum, n_contour, ...
    plotn_chan, plotn_cond, plotn_g, plotn_f, plotn_t, chosen_topochan, ...
    chosen_coh_limit, chosen_coh_limit_diff, f_indiv_hz, chosen_p, plotn_p);

end