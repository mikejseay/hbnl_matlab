function outcellstr=make_timelabels(t_start_ms,t_end_ms)

n_timewins=length(t_start_ms);

outcellstr=cell(n_timewins,1);

for win=1:n_timewins
    outcellstr{win}=sprintf('%d-%d ms',t_start_ms(win),t_end_ms(win));
end

end