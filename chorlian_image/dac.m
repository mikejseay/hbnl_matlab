function [Xcorr, cm] = deautocorrelate(Xdist)
	[nr, nc] = size(Xdist);
	limit = nr - 1;
	Xcorr = zeros(size(Xdist));
	cm = zeros(nr, 1);
	for m = 1:(limit + 1)
		cm(m) = mean(diag(Xdist, m - 1));
		Xcorr = Xcorr + diag(diag(Xdist, m - 1) - cm(m), m - 1);
		if m ~= 1
			Xcorr = Xcorr + diag(diag(Xdist, m - 1) - cm(m), -m + 1);
		end
	end		
