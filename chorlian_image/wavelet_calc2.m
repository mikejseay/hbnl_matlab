function Y = wavelet_calc2(X, scale_vec, srate, cycle_edges)
	% X is a matrix of data (n_samps x n_chans)
	% scale_vec is a vector of scales
    
    if nargin<3
        srate=256;
    end
	
	[n_samps, n_chans] = size(X);
	scale_len = length(scale_vec);
    
    freqs = 2*srate./scale_vec;
    cycles = logspace(log10(cycle_edges(2)),log10(cycle_edges(1)),scale_len); % ./ (2*pi*freqs);
    
	len = pow2(ceil(log2(2 * max(scale_vec) + n_samps)) + 1);
	Xfft = fft(X, len);
	Y = zeros(n_samps, n_chans, scale_len);
	for m = 1:scale_len
		wav_matrix = fast_cwt_wavlet(scale_vec(m), len, n_chans, freqs(m), cycles(m), srate);
		Y(:, :, m) = fast_cwt_wav_mat_calc(Xfft, wav_matrix, ...
			len, n_samps, scale_vec(m));
	end
	
	
function wav_matrix = fast_cwt_wavlet(scale, len, n_chans, freq, cycles, srate)
	% len should be the power of 2 >= 2 * length of signal
	wav = zeros(len, 1);
	wav(:, 1) = mwav(scale, len, cycles, freq, srate);
	wav_matrix = repmat(wav, [1 n_chans]);
		
function w = mwav(scale, len, cycles, freq, srate)
	N = 8 * srate;
	wo = 2 * pi;
	t = linspace(-4, 4, N);
    taper_width=round(scale/2*cycles);
    c = taper_width / N;
    A = 1/sqrt( (cycles/2*pi*freq) *sqrt(pi)); % ?????
	ww = A * exp( -t.^2 ./ (2*c^2)) .* exp(1i * wo * freq .* t);
	www = zeros(len, 1);
	s = 1 + (len - N)/2;
	www(s:s + N - 1) = ww;
	w = fft(www, len);

function Y = fast_cwt_wav_mat_calc(x_fft, wav_matrix, ...
	len, lx, trim)
	Wx = x_fft .* wav_matrix;
	WWx = ifft(Wx, len);
	Y = WWx((1 + len/2):(lx + len/2), :);
	Y([1:trim, (end - trim):end], :) = 0;
