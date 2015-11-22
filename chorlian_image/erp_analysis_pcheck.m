function [Y,name] = erp_analysis_pcheck(data_file, data_file_type, info_file, ...
	param_struct, option, vis)
% stages
% 0) args = data_file, info_file, parameters: wavelet scale, regions
% 1) load file
% get data and trial info
% info file not needed for hdf_binary type
% 2) filter
% 
% 3) do wavelet analysis 
% input = n_samps X n_chans X n_trials
% output = n_samps X n_chans X n_trials x n_freqs for input to state_space
% output = n_samps X n_chans X n_freqs X n_cases for total and evoked
% option = do bootstrapping 
% 
% 4) do state space analysis
% input = n_samps X n_chans X n_trials x n_freqs 
% output = n_segs X n_regions X n_freqs X n_cases
% 
% do bootstrapping 
% how to test for within time_series differences 
% options: 1 wavelet decomp, 2: wavelet coh, 3: state-space
% 4 wavelet cont 
n_default_options = 4;
default_option = true(n_default_options, 1);

%get the subject id
[~,name,~]=fileparts(data_file);

if vis
%set the output path for datachecking figures
check_path=param_struct.pcheck_path;
%check for existence of subdirectory
if ~exist([check_path,'/',name,'/'],'dir')
    mkdir([check_path,'/',name,'/'])
end
subj_ckpath=[check_path,'/',name,'/'];
clear check_path
end

is_mc_bin = false;
if nargin < 5 || isempty(option) || ~islogical(option)
	option = default_option;
else 
	n_options = numel(option);
	if n_default_options < n_options
		option = option(1:n_default_options);
	elseif n_default_options > n_options
		option = [option(:); false(n_default_options - n_options, 1)];
	end
end
clear default_option n_default_options n_options

% stage 1
	Y = [];
	% get parameter from param_struct

	if exist(data_file, 'file') ~= 2
		fprintf(2, 'error: %s does not exist\n', data_file);
		return
	end	

	if strcmpi(data_file_type, 'mc_binary')
		n_chans = param_struct.n_chans;
		n_samps = param_struct.n_samps;
		n_trials = param_struct.n_trials;		
		[data, count] = read_bin(data_file, n_samps * n_trials, n_chans);
		if count == 0
			fprintf(2, 'error reading %s\n', data_file);
		end
		if count < n_samps * n_trials * n_chans
			n_trials_read = floor(count/(n_samps * n_chans));
			fprintf(2, '%s: only read %d trials\n', data_file, n_trials_read);
		end
		is_mc_bin = true;
		if isfield(param_struct, 'reref_chan') && param_struct.reref_chan ~= 0
			dataX = bsxfun(@minus, data, data(:, param_struct.reref_chan));
			dataX(:, param_struct.reref_chan) = ...
				-data(:, param_struct.reref_chan);
			data = dataX; 
			clear dataX;
		end

		trial_info = [];
	
		% if trial_info provided first col is case_num, second is reaction time
		if isnumeric(info_file)  
			trial_info = info_file;
		end
		if exist(info_file, 'file')
			trial_info = load(info_file);
		end
		if isempty(trial_info)
			fprintf(2, '%s: error: no trial info provided\n', data_file);
			return
		end
		nr = size(trial_info, 1);
		if nr < n_trials % zero pad
			fprintf(2, ...
				'%s: error: trial info inconsistent with number of trials ', ...
				data_file);
			fprintf(2, 'zero padding trial info\n');
			nc = size(trial_info, 2);
			new_trial_info = zeros(n_trials, nc);
			new_trial_info(1:nr, :) = trial_info;
			trial_info = new_trial_info;
			clear new_trial_info
		elseif nr > n_trials
			fprintf(2, ...
				'%s: error: trial info inconsistent with number of trials ', ...
				data_file);
			fprintf(2, 'truncating trial info\n');
			trial_info = trial_info(1:n_trials, :);
		end
		trial_type = trial_info(:, param_struct.trial_col);
	elseif strcmpi(data_file_type, 'hdf_binary')
		h1_struct = read_hdf1_dataV7(data_file);
		if isempty(h1_struct) 
			fprintf(2, 'error reading %s\n', data_file);
			return;
		end
		X = h1_struct.data_struct.hdf1_cnt_data;
		if isempty(X)
			fprintf(2, 'no data in %s\n', data_file);
			return;
		else
			data = X';
			clear X
			% at this point data is n_total_samps X n_chans
%CHECK 1 HERE
%return channel standard deviations, time series of mean across channels,
% +/- 2 SD, and "envelope" (greatest excursion among channels)
% can use trimOutlier from eeglab
		end
		% unload variables
		trial_type = check_hdf_trials(h1_struct);
	else
		fprintf(2, 'only mc_binary files can be used at this point\n');
		return;
	end
	
	if isfield(param_struct, 'chan_vec')
		data = data(:, param_struct.chan_vec);
	end
	n_chans = size(data, 2);
    
    % clip 19 points off the end to avoid filtering edge effects
    %here we check to see if ending samples have standard deviation across
    %channels greater than 10000, but it always come out as the terminal 19
    %samples or so (artifact from turning off amplifier most likely)
    data = clip_end_rawdata(data,10000);
    raw_size=size(data);
    
    
    if param_struct.cleanset
        % using ASR
        %put the data into an EEGLAB dataset structure
        EEG=import_hdf_eeglab_inline(data_file,data,h1_struct);
        clear data
        
        if vis
            %visualize with trimoutlier
            [EEG,h] = pop_trimOutlierPlot(EEG);
            export_fig([subj_ckpath,name,'_1raw'],'-png');
            close(h)
        end
        
        % high-pass filter at 1 Hz and visualize again (not used
        EEG2 = pop_eegfiltnew(EEG, [], 1, 1650, true, [], 0);
        
        if vis
            %visualize with trimoutlier
            [EEG2,h] = pop_trimOutlierPlot(EEG2);
            export_fig([subj_ckpath,name,'_2hpf'],'-png');
            close(h)
        end
        
        %save channel locations
        chan_locs=EEG.chanlocs;
        
        % filter and clean using ASR
        EEG = cleanset_eeglab(EEG);
        
        %interpolate any missing channels
        EEG = pop_interp(EEG, chan_locs, 'spherical');
        
        if vis
            %visualize with trimoutlier
            [EEG,h] = pop_trimOutlierPlot(EEG);
            export_fig([subj_ckpath,name,'_3asr'],'-png');
            close(h)
        end
        
        dataF = EEG.data';
        clear EEG
    else
        % filter data without regard to trial structure
        dataF = filter_data(data, ...
            param_struct.lower, param_struct.upper, [], param_struct.rate);
        clear data
    end
    
    %check to make sure there are the same number of data points as before
    if any(size(dataF) ~= raw_size)
        fprintf(2, 'error: %s caused a problem with filtering and/or ASR\n', data_file);
        return
    end
    clear raw_size
    
    % reference the data
    %dataR = reference_data(dataF,'none');
    %clear dataF
    
    % CSD transform the data if applicable
    if isfield(param_struct,'csd_G') && isfield(param_struct,'csd_H')
        dataF = dataF'; % transpose to n_chans x n_samps
        dataR = CSD(dataF,param_struct.csd_G,param_struct.csd_H);
        clear dataF
        dataR = dataR'; % transpose to n_samps x n_chans
    end
    
    %put the data into an EEGLAB dataset structure
    EEG=import_hdf_eeglab_inline(data_file,dataR,h1_struct);
    
    if vis
        %visualize with trimoutlier
        [EEG,h] = pop_trimOutlierPlot(EEG);
        export_fig([subj_ckpath,name,'_4csd'],'-png');
        clear EEG
        close(h)
    end
    
	% do wavelet calcuation
	%n_scales = numel(param_struct.scale);
	%dataW = wavelet_calc(dataR, param_struct.scale);
    % vars to save at this point
    % data_file, data_file_type, is_mc_bin
    % param_struct, option, h1_struct
    % dataW, dataR
    % n_samps, n_chans, n_scales, n_cases, trial_type
    % midpoint_struct=v2struct(data_file, data_file_type, is_mc_bin, param_struct, ...
    % option, h1_struct, dataW, dataR, n_chans, n_scales, trial_type);
    

    %% WAVELET DATA WOULD BE SAVED AT THIS POINT (BUT TOO BIG)
    
	if is_mc_bin
		dataW = reshape(dataW, [n_samps, n_trials, n_chans, n_scales]);
	else
		%[~, n_samps, dataW] = hdf_reshapeW(dataW, h1_struct, param_struct.rate, param_struct.prestim_ms);
		[~, n_samps, dataR] = hdf_reshape(dataR, h1_struct, param_struct.rate, param_struct.prestim_ms); 
		if isempty(dataR)
			fprintf(2, 'invalid time data in %s\n', data_file);
			return
		end
		n_trials = size(dataR, 2);
		if n_trials < numel(trial_type)
			trial_type = trial_type(1:n_trials);
		end
	end 
	% dataR is n_samps X n_trials X n_chans
	if isfield(param_struct, 'paired_trials') && ...
			param_struct.paired_trials > 0
		if param_struct.paired_trials == 2 || param_struct.paired_trials == 4
			if mod(n_trials, 2) == 1
				n_trials = n_trials - 1;
				dataW = dataW(:, 1:(end - 1), :, :);
				dataR = dataR(:, 1:(end - 1), :);
			end
			n_trials = n_trials/2;
			n_samps = n_samps * 2;
			trial_type = [trial_type(1:2:end), trial_type(2:2:end)];
			[bad_trial_rows, ~] = find(trial_type == 0);
			if ~isempty(bad_trial_rows)
				good_trial_rows = setdiff(1:n_trials, bad_trial_rows);
			else
				good_trial_rows = 1:n_trials;
			end
			trial_pair_type = zeros(n_trials, 1);
			[trial_types, ~, trial_pair_type(good_trial_rows)] = ...
				unique(trial_type(good_trial_rows, :), 'rows');
			trial_type = trial_pair_type;
			dataW = reshape(dataW, [n_samps, n_trials, n_chans, n_scales]);
			dataR = reshape(dataR, [n_samps, n_trials, n_chans]);
		end
		if param_struct.paired_trials == 4
			trial_pairs = [trial_type(1:(end - 1)), trial_type(2:end)];
			[bad_trial_rows, ~] = find(trial_pairs == 0);
			if ~isempty(bad_trial_rows)
				good_trial_rows = setdiff(1:(n_trials - 1), bad_trial_rows);
			else
				good_trial_rows = 1:(n_trials - 1);
			end			
			trial_pair_type = zeros(n_trials - 1, 1);
			[trial_types, ~, trial_pair_type(good_trial_rows)] = ...
				unique(trial_pairs(good_trial_rows, :), 'rows');
			trial_type(2:end) = trial_pair_type;
			trial_type(1) = 0;
		end
	else
		trial_types = trial_type;
    end
    if isfield(param_struct, 'trial_init_type')
        if param_struct.trial_init_type==90
            etable=h1_getbehav(h1_struct,param_struct.trial_init_type);
            [param_struct.trial_mat,param_struct.case_label,respRT] = ...
                behavinds_sog(etable);
        end
    end
	n_cases = max(trial_type); %is this advisable?
	if ~isfield(param_struct, 'trial_mat')
		trial_mat = false(n_trials, n_cases);
		for m = 1:n_cases
			trial_mat(:, m) = trial_type == m;
		end
		if isfield(param_struct, 'case_vec') %in other words, if there is a
            % case vector specified, take only those, otherwise take all
            % cases
			trial_mat = trial_mat(:, param_struct.case_vec);
		end
		zero_trials = find(sum(trial_mat) == 0);
		if ~isempty(zero_trials)
			trial_mat(:, zero_trials) = [];
		end
	else
		trial_mat = param_struct.trial_mat;
	end
	n_cases = size(trial_mat, 2);
	%check to see if last trial was removed and reflect that change in the
	%trial_mat
    %if length(trial_indexW) ~= size(trial_mat,1)
    %    trial_mat = trial_mat(1:(end - 1), :);
    %end
    
% end of stage 1	
% filter and calculate CWT

	% eliminate first and last trials
	
	%dataW = dataW(:, 2:(end - 1), :, :); %takes forever???
	dataR = dataR(:, 2:(end - 1), :);
    trial_mat = trial_mat(2:(end - 1), :);
    if isfield(param_struct, 'artf_chan_vec')
        artf_chan_vec = param_struct.artf_chan_vec;
    else
        artf_chan_vec = 1:n_chans;
    end
    if isfield(param_struct, 'artf_samp_vec')
        artf_samp_vec = param_struct.artf_samp_vec;
    elseif isfield(param_struct, 'artf_ms_vec')
        artf_samp_vec = round(param_struct.artf_ms_vec/1000*...
            param_struct.rate);
        artf_samp_vec=artf_samp_vec(1):artf_samp_vec(2);
    else
        if isfield(param_struct, 'prestim_ms')
            artf_samp_vec = (49:166) + ...
               round(param_struct.prestim_ms/1000 * param_struct.rate);
        else
            artf_samp_vec = 49:166; %????
        end
    end
    %artifact rejection by trial
    %if CSD-transformed
    if isfield(param_struct,'csd_G')
    trial_mat_orig=trial_mat;
    for thresh=150:50:1500
        trial_mat = check_artifact(dataR(artf_samp_vec, :, :), trial_mat_orig, ...
            artf_chan_vec, thresh); %artifact rejection
        prop_kept=sum(sum(trial_mat))/sum(sum(trial_mat_orig));
        if prop_kept > .75
            fprintf(2, 'threshold used was %d\n', thresh);
            break
        end
    end
    %if not CSD-transformed but has a thresh (usually 150)
    elseif isfield(param_struct, 'thresh')
        thresh=param_struct.thresh;
        trial_mat = check_artifact(dataR(artf_samp_vec, :, :), trial_mat, ...
            artf_chan_vec, thresh); %artifact rejection
    %if not given enough info about how to reject?    
    else
        thresh=500;
        trial_mat = check_artifact(dataR(artf_samp_vec, :, :), trial_mat, ...
            artf_chan_vec, thresh); %artifact rejection
    end
	if isfield(param_struct, 'max_trials') && ...
			isfield(param_struct, 'min_trials')
		trial_mat = fix_trial_mat(trial_mat, param_struct.min_trials, ...
			param_struct.max_trials);
	end
	%wavelet_evk = zeros(n_samps, n_chans, n_scales, n_cases);
	%wavelet_tot = zeros(n_samps, n_chans, n_scales, n_cases);
	mean_data = zeros(n_samps, n_chans, n_cases);
    % check to see if there are enough trials in the 
	for m = 1:n_cases
		% keep phase info
		%wavelet_evk(:, :, :, m) = ...
		%	squeeze(mean(dataW(:, trial_mat(:, m), :, :), 2));
		%wavelet_tot(:, :, :, m) = ...
		%	squeeze(mean(abs(dataW(:, trial_mat(:, m), :, :)), 2));
		mean_data(:, :, m) = squeeze(mean(dataR(:, trial_mat(:, m), :), 2));
	end
	%Y.wavelet_evk = wavelet_evk;
	%Y.wavelet_tot = wavelet_tot;
    Y.trial_data = dataR;
	Y.mean_data = mean_data;
	Y.trial_mat = trial_mat;
	Y.n_trials = sum(trial_mat);
    Y.etable = etable;
    if isfield(param_struct,'getbehav')
        if param_struct.getbehav && param_struct.trial_init_type==90
            
        [avgbet_prevoutcome,stdbet_prevoutcome,crit,crit_prevoutcome,avgbet,resp_mat, ...
 prevoutcome_mat,avgbet2_prevoutcome,stdbet2_prevoutcome] = behav_sog(etable);

        behav_data=v2struct(avgbet_prevoutcome,stdbet_prevoutcome,crit,crit_prevoutcome, ...
            avgbet,resp_mat, prevoutcome_mat,avgbet2_prevoutcome,stdbet2_prevoutcome,respRT);
        
        Y.behav_data = behav_data;
        
        end
    end
	
	if sum(option) == 0
% 		Y.dataW = dataW;
		
		return		
	end
	if isfield(param_struct, 'n_perms')
		n_perms = param_struct.n_perms;
		seed = param_struct.seed;
	else
		n_perms = 0;
		seed = 0;
	end
% % do wavelet coherence analysis
	if option(1)
		coherence_pairs = param_struct.coherence_pairs;
		if isfield(param_struct, 'chan_vec') % check for consistency with pairs
			[U, ~, Un] = unique(param_struct.coherence_pairs);
			if isequal(U, param_struct.chan_vec)
				coherence_pairs = reshape(Un, size(param_struct.coherence_pairs));
			elseif all(ismember(U, 1:n_chans))
				coherence_pairs = param_struct.coherence_pairs;
			else
				fprintf(2, ...
					'Inconsistent specification of chan_vec and coherence_pairs\n');
				coherence_pairs = [];
			end
		end
		if ~isempty(coherence_pairs)
			n_coh_pairs = size(param_struct.coherence_pairs, 1);
			if n_perms > 0
				coh_results{1} = zeros(n_samps, n_scales, n_cases, n_coh_pairs);
				coh_results{2} = zeros(n_samps, n_scales, n_cases, n_coh_pairs);
				for m = 1:n_coh_pairs
					[coh_results{1}(:, :, :, m),  coh_results{2}(:, :, :, m)] = ...
						pairwise_coherence_calc(dataW(:, :, coherence_pairs(m, :), :), ...
							trial_mat, n_perms, seed);
				end
			else
				coh_results = zeros(n_samps, n_scales, n_cases, n_coh_pairs);
				for m = 1:n_coh_pairs
					coh_results(:, :, :, m) = ...
						pairwise_coherence_calc(dataW(:, :, coherence_pairs(m, :), :), ...
							trial_mat, n_perms, seed);
				end		
			end
		else
			coh_results = [];
		end
		Y.coh_results = coh_results;
	end

% % do wavelet cont coherence analysis
	if option(2)
		coherence_pairs = param_struct.coherence_pairs;
		if isfield(param_struct, 'chan_vec') % check for consistency with pairs
			[U, ~, Un] = unique(param_struct.coherence_pairs);
			if isequal(U, param_struct.chan_vec)
				coherence_pairs = reshape(Un, size(param_struct.coherence_pairs));
			elseif all(ismember(U, 1:n_chans))
				coherence_pairs = param_struct.coherence_pairs;
			else
				fprintf(2, ...
					'Inconsistent specification of chan_vec and coherence_pairs\n');
				coherence_pairs = [];
			end
		end
		n_coh_pairs = size(param_struct.coherence_pairs, 1);
		if n_perms > 0
			for m = 1:n_coh_pairs
				[U1, U2] = ...
					pairwise_cont_coh_calc(dataW(:, :, coherence_pairs(m, :), :), ...
						trial_mat, 32, 16, n_perms, seed);
				if m == 1
					Usize = size(U1);
					cont_coh_results = zeros([Usize, n_coh_pairs, 2]);
				end
				cont_coh_results(:, :, :, m, 1) = U1;
				cont_coh_results(:, :, :, m, 2) = U2;
			end
		else
			for m = 1:n_coh_pairs
				U = ...
					pairwise_cont_coh_calc(dataW(:, :, coherence_pairs(m, :), :), ...
						trial_mat, 32, 16, n_perms, seed);
				if m == 1
					Usize = size(U);
					cont_coh_results = zeros([Usize, n_coh_pairs]);
				end			
			cont_coh_results(:, :, :, m) = U;
			end
		end
		Y.cont_coh_results = cont_coh_results;
	end	
	
% % do state-space analysis
	if option(3) && isfield(param_struct, 'regions')
		n_regions = numel(param_struct.regions);
		n_times = param_struct.n_times;
		overlap = param_struct.overlap;
		regional_sync_stat = cell(n_regions, n_cases, n_scales);
		dataWW = permute(dataW, [1, 3, 2, 4]);
		regional_sync = [];
		for m = 1:n_regions
			rg = param_struct.regions{m};
			for n = 1:n_cases
				for r = 1:n_scales
					[U, ~, regional_sync_stat{m, n, r}] = ...
						svd_entropy(dataWW(:, rg, trial_mat(:, n), r), ...
						n_times, overlap, n_perms, seed);
					if isempty(regional_sync)
						nr = size(U, 1);
						regional_sync = zeros(nr, n_scales, n_cases, n_regions);
					end
					regional_sync(:, r, n, m) = U;
				end
			end
		end
		clear dataWW
		Y.regional_sync = regional_sync;
		Y.regional_sync_stat = regional_sync_stat;
    end
    
% % do phase-amp analysis
	if option(4) && isfield(param_struct, 'phase_amp_chan_vec')
		chan_vec = param_struct.phase_amp_chan_vec;
		phase_amp = cell(n_cases, 1);
		for m = 1:n_cases
			phase_amp{m} = calc_phase_amp(dataW(:, :, trial_mat(:, n), :), ...
				chan_vec, scale_mat, range_mat);
		end
		Y.phase_amp = phase_amp;
	end

