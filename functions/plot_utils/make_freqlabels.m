function outcellstr=make_freqlabels(f_start_hz,f_end_hz)

n_freqwins=length(f_start_hz);

outcellstr=cell(n_freqwins,1);

for win=1:n_freqwins
    outcellstr{win}=sprintf('%1.1f-%1.1f hz',f_start_hz(win),f_end_hz(win));
end

end