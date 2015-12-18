function [Y, filt_coefs] = filter_data(X, a, b, filt_coefs, rate)
	[n_samps, n_chans] = size(X);
	if (nargin == 3) || isempty(filt_coefs)
		if nargin == 5
			filt_coefs = make_filt_coefs(a, b, rate);
		else
			filt_coefs = make_filt_coefs(a, b);
		end
	end
	[B, nfft, L, nx] = fftfilt_start(filt_coefs, X);
	Y = fftfilt_finish(B, X, nfft, L, nx);
	XXX = fftfilt_finish(B, flipud(Y),  nfft, L, nx);
	Y = flipud(XXX);

function filt_coefs = make_filt_coefs(a, b, rate)
	if nargin == 2
		rate = 256;
	end
	c = (a + b)/2; w = (b - a)/2;
	f = [0. (c - w)/(rate/2), (c - w)/(rate/2), (c + w)/(rate/2), ...
		(c + w)/(rate/2), 1];
	m = [0        0        1            1           0      0];
	filt_coefs = firls(rate - 1, f, m);	
