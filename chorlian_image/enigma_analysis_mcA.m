function output = enigma_analysis_mc(file, chan_idx)
	if nargin == 1 || isempty(chan_idx)
		chan_idx = [19, 1, 15];
	end
	output = [];
	[Y, trial_idx] = enigma_pre_analysis(file, chan_idx);
	if isempty(Y) || sum(trial_idx < 128)
		return;
	end
	cFreq = [.5, 3.5; 4, 7.5; 8, 12.5; 13, 30; .5, 30];
	peak_range = [7, 14];
	
	cWinSec = 2;      % use 2 second length windows
	cRate = 256;      % 250 was our sampling rate
	cWinSize = cRate * cWinSec; % FFT resolution will be .5 Hz
	cNChans    = size(Y, 2);  % set these correctly
	cNsampsWin = cWinSize * cRate;
	
	fs = linspace(0, cRate/2, cWinSize/2 + 1);
	fs2 = linspace(0,cRate/2, cWinSize + 1);
	hanning_window = repmat(hanning(cWinSize), 1, cNChans);
	n_trials = numel(trial_idx);
	trial_start = 1:256:2^16;
	trial_finish = trial_start + 255;
	trial = 1;
	windows_used = 0;
	P = zeros(cWinSize, numel(chan_idx));
	while trial < n_trials - 1
		if trial_idx(trial) && trial_idx(trial + 1)
			F = fft(hanning_window .* ...
				Y(trial_start(trial):trial_finish(trial + 1), :));
			P = P + (F .* conj(F))/cNsampsWin;
			windows_used = windows_used + 1;
		end
		trial = trial + 1;
	end
	if windows_used < 128
		return;
	end
	PW = P(1:(end/2) + 1, :) / windows_used;
	PW(2:(end - 1), :) = 2 * PW(2:(end - 1), :);
	PowF = zeros(5, 3);
	for m = 1:5
		f_idx = fs >= cFreq(m, 1) & fs <=  cFreq(m, 2);
		PowF(m, :) = mean(log(PW(f_idx, :)));
	end
	output.PowF = PowF;
% 	f_idx = fs2 >= peak_range(1) & fs <= peak_range(2);	
