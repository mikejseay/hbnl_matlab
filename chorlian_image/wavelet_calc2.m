function Y = wavelet_calc2(X, scale_vec, cycles_vec)
	% X is a matrix of data
	% scale_vec is a vector of scales
    % cycles_vec is a corresponding vector of cycles
	% procedure: fft of X, create fft of wave lets
	[n_samps, n_chans] = size(X);
% 	scale_vec = flipud(sort(scale(:)));
	scale_len = length(scale_vec);
	len = pow2(ceil(log2(2 * max(scale_vec) + n_samps)) + 1);
	Xfft = fft(X, len);
	Y = zeros(n_samps, n_chans, scale_len);
	for m = 1:scale_len
		wav_matrix = fast_cwt_wavlet(scale_vec(m), cycles_vec(m), len, n_chans);
		Y(:, :, m) = fast_cwt_wav_mat_calc(Xfft, wav_matrix, ...
			len, n_samps, scale_vec(m));
	end
	
	
function wav_matrix = fast_cwt_wavlet(scale, cycle, len, n_chans)
	% len should be the power of 2 >= 2 * length of signal
	wav = zeros(len, 1);
	wav(:, 1) = mwav(scale, cycle, len);
	wav_matrix = repmat(wav, [1 n_chans]);
		
function w = mwav(scale, cycle, len)
    width = round( 2 .* cycle );
	N = width .* scale; %length in samples
	wo = 2 * pi;
	t = linspace(-width, width, N); %analytical time
    cycle_const = (cycle+1) ./ (2 * pi);
    gauss = exp(-(t .^ 2)./(2.*cycle_const.^2));
    ww = sqrt(1/(scale.*cycle_const)) .* exp(1i * wo .* t) .* gauss;
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
