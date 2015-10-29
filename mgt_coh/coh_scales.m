if ~exist('tfparams','var')

    %time scaling (in pts for now, can convert from a ms declaration in
    %parameters)
    t_start=1;
    t_end=imp.maxtimepts;
    %define absolute ms edges
    ms_start=-param_struct.prestim_ms;
    ms_end=ms_start + ( imp.maxtimepts * (1000 / param_struct.timerate) );
    t_ms=linspace(ms_start, ms_end, imp.maxtimepts+1);
    %clip surplus points beyond declared end point
    t_ms(imp.maxtimepts+1:end)=[];
    ms_tickint=200;
    t_xtick_ms=(ms_start+mod(ms_start,ms_tickint)): ...
        ms_tickint:(ms_end-mod(ms_end,ms_tickint));
    t_xtick=zeros(length(t_xtick_ms),1);
    for tick=1:length(t_xtick_ms)
        [~,i]=min(abs(t_ms-t_xtick_ms(tick)));
        t_xtick(tick)=i;
    end
    [~,t_zero]=min(abs(t_ms));

    %freqs
    freqs=param_struct.rate*2./param_struct.scale;
    f_ytick=[20 17 15 13 8 4];
    f_label=round(freqs(f_ytick));
    f_ytick=imp.maxfreqs-f_ytick+1;
    f_nlabel=cellstr(int2str(round(freqs)'));
    f_color=distinguishable_colors(length(freqs));

else
    
    %times
    t_start=1; t_end=tfparams.n_times;
    t_ms=tfparams.times;
    ms_tickint=200;
    t_xtick_ms=(t_ms(1)+mod(t_ms(1),ms_tickint)): ...
        ms_tickint:(t_ms(end)-mod(t_ms(end),ms_tickint));
    t_xtick=zeros(length(t_xtick_ms),1);
    for tick=1:length(t_xtick_ms)
        [~,i]=min(abs(t_ms-t_xtick_ms(tick)));
        t_xtick(tick)=i;
    end
    [~,t_zero]=min(abs(t_ms));
    
    %freqs
    freqs=flip(tfparams.freqs);
    f_ytick=[28 26 24 17 7];
    f_label=round(freqs(f_ytick));
    f_ytick=imp.maxfreqs-f_ytick+1;
    f_nlabel=cellstr(int2str(round(freqs)'));
    f_color=distinguishable_colors(length(freqs));
    
end

%conditions - OK
cond_label=param_struct.case_label;

%subjects - OK
s_label=cellstr(int2str((1:imp.s_valid)'));
s_color=distinguishable_colors(imp.s_valid);

%channel pairs - OK
p_label=cell(imp.maxpairs,1);
for pair=1:imp.maxpairs
    p_label{pair}=[chan_locs(param_struct.coherence_pairs(pair,1)).labels,'-',chan_locs(param_struct.coherence_pairs(pair,2)).labels];
end
p_color=distinguishable_colors(length(p_label));

%channels -OK
chan_label=cell(length(chan_locs),1);
for chan=1:imp.maxchans
    chan_label{chan}=chan_locs(chan).labels;
end
chan_color=distinguishable_colors(length(chan_label));

%phi (for making figs)
phi = GoldenRatio();

%pack into scale struct
scl=v2struct(chan_color, chan_label, ...
    cond_label, ...
    f_color, f_label, f_nlabel, f_ytick, freqs, ...
    p_color, p_label, ...
    s_color, s_label, ...
    t_end, t_start, t_ms, t_xtick, t_xtick_ms, t_zero,...
    phi);

clear chan_color chan_label cond_label f_color f_label f_nlabel f_ytick ...
    freqs p_color p_label s_color s_label t_end t_start t_ms t_xtick t_xtick_ms t_zero

clear ms_start ms_end ms_tickint tick c i cond_inds chan chosen_s cond ...
    cond_label_seed freq pair phi