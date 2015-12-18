function Y = wavelet_calc(X, scale_vec)
	% X is a matrix of data
	% scale is a vector of scales
	% procedure: fft of X, create fft of wave lets
	[n_samps, n_chans] = size(X);
% 	scale_vec = flipud(sort(scale(:)));
	scale_len = length(scale_vec);
	len = pow2(ceil(log2(2 * max(scale_vec) + n_samps)) + 1);
	Xfft = fft(X, len);
	Y = zeros(n_samps, n_chans, scale_len);
	for m = 1:scale_len
		wav_matrix = fast_cwt_wavlet(scale_vec(m), len, n_chans);
		Y(:, :, m) = fast_cwt_wav_mat_calc(Xfft, wav_matrix, ...
			len, n_samps, scale_vec(m));
	end
	
	
function wav_matrix = fast_cwt_wavlet(scale, len, n_chans)
	% len should be the power of 2 >= 2 * length of signal
	wav = zeros(len, 1);
	wav(:, 1) = mwav(scale, len);
	wav_matrix = repmat(wav, [1 n_chans]);
		
function w = mwav(scale, len)
	N = 4 * scale; %length in samples
	wo = 2 * pi;
	t = linspace(-4, 4, N); %analytical time
	ww = sqrt(1/scale) * exp(1i * wo .* t) .* exp(-(t .^ 2)/2); %gives roughly 6 cycles
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
