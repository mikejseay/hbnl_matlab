function [Y, trial_idx] = enigma_artf_check(X,diffthresh)	
	[ns, nc] = size(X);
	Y = detrend(X);
	nt = floor(ns/256);
	YY = reshape(Y(1:(nt * 256), :), [256, nt, nc]);
	trial_idx = true(nt, 1);
	% check for artifact 
	for m = 1:nc
		for n = 1:nt
			if trial_idx(n)
				XX = detrend(YY(:, n, m));
				stdX = std(XX);
				diffX = max(XX) - min(XX);
				if stdX < 1 || diffX > diffthresh 
					trial_idx(n) = false;
				end
			end
		end
	end
