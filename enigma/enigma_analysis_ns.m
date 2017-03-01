function output = enigma_analysis_ns(file, chan_idx)
	if nargin == 1 || isempty(chan_idx)
		chan_idx = [16, 31, 30];
	end
	h1_struct = read_hdf1_dataV7(file); 
	X = h1_struct.data_struct.hdf1_cnt_data;
	data = X(chan_idx, :)';
	n_samps = size(data, 1);
	clear X
	diff_n_samps = n_samps - (254 * 256);
	if diff_n_samps > 0
		start = floor(diff_n_samps/2);
		finish = start + (254 * 256) - 1;
		data = data(start:finish, :);
	end 
	
	dataF = filter_data(detrend(data), .2, 40, [], 256);
	[Y, trial_idx] = enigma_artf_check(dataF,125);
	clear dataF
	if isempty(Y) || sum(trial_idx) < 128
		output=sum(trial_idx);
        %subplot(1,2,2); plot(data)
    else
	output = enigma_analysis(Y, trial_idx);
    end
end
