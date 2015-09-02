function output = enigma_analysis_mc(file, chan_idx)
	if nargin == 1 || isempty(chan_idx)
		chan_idx = [20, 1, 16];
	end
	output = [];
	X = read_bin(file, 256 * 256, 21);
   if size(X, 2) ~= 21
      return
   end
	data = X(257:(end - 256), chan_idx);
	clear X
	dataF = filter_data(detrend(data), .2, 40, [], 256);
	[Y, trial_idx] = enigma_artf_check(dataF,125);
	clear dataF
	if isempty(Y) || sum(trial_idx) < 128
		output=sum(trial_idx);
        %figure; subplot(1,2,1); plot(data)
    else
	output = enigma_analysis(Y, trial_idx);
    end
end
