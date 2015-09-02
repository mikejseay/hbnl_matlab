function [Y, trial_idx] = enigma_pre_analysis(file, chan_idx)
	% chan_idx for enigma = [19, 1, 15]; = CZ, O1, O2
	if nargin == 1 || isempty(chan_idx)
		chan_idx = [19, 1, 15];
	end
	data = read_bin(file, 256 * 256, 20);
   if size(data, 2) ~= 20
      Y = []; trial_idx =[];
      return
   end
	nc = numel(chan_idx);
	Y = data(:, chan_idx);
	YY = reshape(Y, [256, 256, nc]); 
	trial_idx = true(256, 1);
	% check for artifact 
	for m = 1:nc
		for n = 1:256
			if trial_idx(n)
				X = YY(:, n, m);
				stdX = std(X);
				diffX = max(X) - min(X);
				if stdX < 1 || diffX > 100 
					trial_idx(n) = false;
				end
			end
		end
	end
