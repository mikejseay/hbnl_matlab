function [scl,chan_locs]=coh_scales(opt, imp, ms_start_plot)

if nargin<3
    ms_start_plot = -500;
end

if length(opt) > 1
    oopt = opt;
    opt = opt(1);
    n_opts = length(oopt);
else
    n_opts = 1;
end

%time scaling 
t_start=1; t_end=imp.tfmaxtimepts;
if isfield(opt, 'prestim_ms')
    ms_start=-opt.prestim_ms;
else
    ms_start=opt.epoch_lims(1);
end
ms_end=ms_start + ( imp.tfmaxtimepts * (1000 / imp.tftimerate) );
t_ms=linspace(ms_start, ms_end, imp.tfmaxtimepts+1);

%clip surplus points beyond declared end point
t_ms(imp.tfmaxtimepts+1:end)=[];

ms_tickint=200;
[~, t_start_ind] = min(abs(t_ms-ms_start_plot));
t_xtick_ms=(ms_start_plot-rem(ms_start_plot,ms_tickint)): ...
    ms_tickint:(ms_end-mod(ms_end,ms_tickint));
t_xtick=zeros(length(t_xtick_ms),1);
for tick=1:length(t_xtick_ms)
    [~,i]=min(abs(t_ms-t_xtick_ms(tick)));
    t_xtick(tick)=i;
end
[~,t_zero]=min(abs(t_ms));

% time scaling for ERPs
t_start_erp=1; t_end_erp=imp.erpmaxtimepts;
%ms_start=-opt.prestim_ms;
ms_end=ms_start + ( imp.erpmaxtimepts * (1000 / imp.erptimerate) );
t_ms_erp=linspace(ms_start, ms_end, imp.erpmaxtimepts+1);

%clip surplus points beyond declared end point
t_ms_erp(imp.erpmaxtimepts+1:end)=[];

[~, t_start_ind_erp] = min(abs(t_ms_erp-ms_start_plot));
t_xtick_erp=zeros(length(t_xtick_ms),1);
for tick=1:length(t_xtick_ms)
    [~,i]=min(abs(t_ms_erp-t_xtick_ms(tick)));
    t_xtick_erp(tick)=i;
end
[~,t_zero_erp]=min(abs(t_ms_erp));

t_xtick_ms = t_xtick_ms ./ 1000;

%freqs
freqs=opt.rate*2./opt.wavelet_scales;

% these are custom tick skins for certain commonly used scales
if length(freqs)==20
    f_ytick=[20 17 15 13 8 4];
elseif length(freqs)==31
    f_ytick=[31 27 22 19 16 11 4 1];
elseif length(freqs)==16
    f_ytick=[16 14 12 10 8 6 2 1];
end

f_label=round(freqs(f_ytick)*10)/10;
f_label2=round(freqs(f_ytick));
f_ytick=imp.maxfreqs-f_ytick+1;
f_nlabel=cellstr(num2str(freqs',2));
f_color=distinguishable_colors(length(freqs));

%conditions
if n_opts > 1
    cond_label = cell(n_opts*length(opt.case_vec), 1);
    cd = 0;
    for o=1:n_opts
        for c=1:length(oopt(o).case_vec)
            cd = cd + 1;
            cond_label{cd} = [oopt(o).outsuffix, '_', oopt(o).case_label{c}];
        end
    end
else
    cond_label=opt.case_label;
end

%subjects
s_label=cellstr(int2str((1:imp.s_valid)'));
s_color=distinguishable_colors(imp.s_valid);

%load locations
load(opt.coords_file)

%channel pairs - OK
p_label=cell(imp.maxpairs,1);
for pair=1:imp.maxpairs
    p_label{pair}=[chan_locs(opt.coherence_pairs(pair,1)).labels,'-',chan_locs(opt.coherence_pairs(pair,2)).labels];
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
    f_color, f_label, f_label2, f_nlabel, f_ytick, freqs, ...
    p_color, p_label, ...
    s_color, s_label, ...
    t_end, t_start, t_ms, t_xtick, t_xtick_ms, t_zero,...
    phi, t_ms_erp, t_xtick_erp, t_zero_erp, ...
    t_start_erp, t_end_erp, ...
    t_start_ind, t_start_ind_erp);

end