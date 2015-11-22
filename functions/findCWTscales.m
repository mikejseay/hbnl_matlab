function scales=findCWTscales(n_samps,srate,freq_lims,padratio)

%srate=256;
%n_samps=384;
%freq_lims=[2 30];
%padratio=2;
%cycles=[3 0.5];

win=max(pow2(nextpow2(n_samps)-3),4);
n_freqs=win/2*padratio+1;
tmpfreqs=linspace(0,srate/2,n_freqs);
tmpfreqs=tmpfreqs(2:end);
n_freqs = length(tmpfreqs( intersect( find(tmpfreqs >= freq_lims(1)), find(tmpfreqs <= freq_lims(2)))));

freqs=linspace(log(freq_lims(1)),log(freq_lims(2)),n_freqs);
freqs=exp(freqs);
freqsb=fliplr(freqs);

scales=round(2*srate./freqsb);

end

%cycles = [ cycles(1) cycles(1)*freqs(end)/freqs(1)*(1-cycles(2))];

%win = dftfilt3(freqs,cycles,srate, 'cycleinc', 'log');