
function [YY, Z, B] = pairwise_coherence_calc_noevents(data, trial_mat, n_perms, seed, type)
	% data is n_samps X n_trials X 2(chans) X n_scales
    
    if nargin<5
        type='pure';
    end
    
	[n_samps, ~, ~, n_scales] = size(data);
	X = squeeze(data(:, :, 1, :) .* conj(data(:, :, 2, :)));
    % X is n_samps X n_trials X n_scales
    
	n_cases = size(trial_mat, 2);
    
	YY = zeros(n_scales, n_cases);
	
	Z = []; B = [];
    if strcmpi(type,'weighted')
        %for n = 1:n_cases
        %    trial_vec = trial_mat(:, n);
        %    YY(:, :, n) = squeeze(sum( X(:, trial_vec, :),1 ) ./ ...
        %        sum(abs(X(:, trial_vec, :)),1 ));
        %end
    elseif strcmpi(type,'pure')
        Xdem = squeeze(data(:, :, 1, :) .* data(:, :, 2, :)); %no conjugate
        Y = squeeze(sum(X, 1) ./ sum(abs(Xdem), 1));
        % Y is now n_trials x n_scales
        
        for n = 1:n_cases
            trial_vec = trial_mat(:, n);
            YY(:, n) = squeeze( mean(Y(trial_vec, :),1) );
        end
    end
    
%{    
	if n_perms > 0 %FIX THIS PART IF YOU WANT, I DIDN'T MAKE IT WORK RIGHT YET
		if seed > 0
			rs = RandStream('mt19937ar');
			RandStream.setGlobalStream(rs);
			reset(rs, seed);
		end	
	
		B = zeros(n_samps, n_scales, n_perms);
		Bmean = zeros(n_samps, n_scales, n_cases);
		Bstd = zeros(n_samps, n_scales, n_cases);
		Z = zeros(n_samps, n_scales, n_cases);
		for n = 1:n_cases
			trial_vec = find(trial_mat(:, n));
			n_trials = numel(trial_vec);
			for m = 1:n_perms
				p = randperm(n_trials);
				XX = squeeze(data(:, trial_vec, 1, :) .* ...
					conj(data(:, trial_vec(p), 2, :)));
				B(:, :, m) = squeeze(abs(sum(XX, 2)) ./ sum(abs(XX), 2));
			end
			B = atanh(B);
			BBmean = mean(B, 3);
			BBstd = std(B, 1, 3);
			Z(:, :, n) = (atanh(abs(Y(:, :, n))) - BBmean) ./ BBstd;
% 			Z(:, :, n) = (abs(Y(:, :, n)) - BBmean) ./ BBstd;
			Bmean(:, :, n) = BBmean;
			Bstd(:, :, n) = BBstd;
		end
		clear B
		B.mean = Bmean;
		B.std = Bstd;
	end
%}

