function y = fftfilt_finish(B, x, nfft, L, nx)	
	y = zeros(size(x));
	istart = 1;
	while istart <= nx
	    iend = min(istart + L - 1, nx);
	    if (iend - istart) == 0
	        X = x(istart(ones(nfft, 1)), :);  % need to fft a scalar
	    else
	        X = fft(x(istart:iend, :), nfft);
	    end
	    Y = ifft(X .* B);
	    yend = min(nx, istart + nfft - 1);
	    y(istart:yend, :) = y(istart:yend, :) + Y(1:(yend - istart + 1), :);
	    istart = istart + L;
	end
	y = real(y);

