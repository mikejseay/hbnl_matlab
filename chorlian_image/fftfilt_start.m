function [B, nfft, L, nx] = fftfilt_start(b, x)
	[m,n] = size(x);
	if m == 1
	    x = x(:);    % turn row into a column
	end
	nx = size(x,1);
	if min(size(b))>1
	   if (size(b,2) ~= size(x,2)) && (size(x,2)>1) && (size(b,2) > 1)
	      error('Filter matrix B must have same number of columns as X.')
	   end
	else
	   b = b(:);   % make input a column
	end
	nb = size(b,1);
	

    if nb >= nx     % take a single FFT in this case
        nfft = 2^nextpow2(nb+nx-1);
        L = nx;
    else
        fftflops = [ 18 59 138 303 660 1441 3150 6875 14952 32373 69762 ...
       149647 319644 680105 1441974 3047619 6422736 13500637 28311786 59244791];
        n = 2.^(1:20); 
        validset = find(n>(nb-1));   % must have nfft > (nb-1)
        n = n(validset); 
        fftflops = fftflops(validset);
        % minimize (number of blocks) * (number of flops per fft)
        L = n - (nb - 1);
        [dum,ind] = min( ceil(nx./L) .* fftflops );
        nfft = n(ind);
        L = L(ind);
    end

	B = fft(b,nfft);
	if length(b)==1,
	     B = B(:);  % make sure fft of B is a column (might be a row if b is scalar)
	end
	if size(b,2)==1
	    B = B(:,ones(1,size(x,2)));  % replicate the column B 
	end



