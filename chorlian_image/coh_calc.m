function Y = coh_calc(data_file, opt)
% pre-processes EEG data, applies wavelet transform, and calculates
% coherence (ITC and ERPCOH).

% get the subject id
[~,name,ext]=fileparts(data_file);

if strcmpi(ext,'.h1')
    data_file_type = 'hdf_binary';
else
    data_file_type = 'mc_binary';
end

if opt.vis
    %set the output path for datachecking figures
    if isfield(opt,'vis_path')
        check_path=opt.vis_path;
    else
        check_path=opt.output_dir;
    end
    %check for existence of subdirectory
    if ~exist([check_path,'/',name,'/'],'dir')
        mkdir([check_path,'/',name,'/'])
    end
    subj_ckpath=[check_path,'/',name,'/'];
    clear check_path
end

% initialize output structure
Y = [];
    
% check for existence of needed data (cleaned and/or CSD-transformed)
% build a suffix based on presence of cleaning and CSD
subdir1 = name(1:3);
check_suffix='';
subdir2 = '';
if isfield(opt,'cleanset')
    if opt.cleanset
        check_suffix=[check_suffix,'_clean'];
        subdir2 = [subdir2, 'clean'];
        if isfield(opt, 'burst') && isempty(opt.burst)
            check_suffix=[check_suffix,'_burst'];
            subdir2 = [subdir2, '_burst'];
        end
    end
end
if isfield(opt, 'eog_rem') && ~isempty(opt.eog_rem)
    check_suffix=[check_suffix, '_', opt.eog_rem];
    subdir2 = [subdir2, '_', opt.eog_rem];
end
if isfield(opt,'csd_G') && ~isempty(opt.csd_G) && isfield(opt,'csd_H') && ~isempty(opt.csd_H)
    check_suffix=[check_suffix,'_CSD'];
    subdir2 = [subdir2, '_CSD'];
end
clean_name = [name, check_suffix, '.mat'];
if isfield(opt,'cleandir') && ~isempty(opt.cleandir)
    cleanfile_dir = strjoin( {opt.cleandir, subdir1, subdir2}, filesep);
    cleanfile_path=fullfile(cleanfile_dir, clean_name);
else
    cleanfile_dir=[];
    cleanfile_path=[];
    fprintf('no directory for cleaned h1 recordings was specified\n');
end

% if the data hasn't already been cleaned
if ~exist(cleanfile_path,'file')
%%    
    
    % if the h1 doesn't exist, error
    if exist(data_file, 'file') ~= 2
        fprintf(1, 'error: %s does not exist\n', data_file);
        return
    end
    % if working with masscomp binary data
    if strcmpi(data_file_type, 'mc_binary')
		n_chans = opt.n_chans;
		n_samps = opt.n_samps;
		n_trials = opt.n_trials;		
		[data, count] = read_bin(data_file, n_samps * n_trials, n_chans);
		if count == 0
			fprintf(1, 'error reading %s\n', data_file);
		end
		if count < n_samps * n_trials * n_chans
			n_trials_read = floor(count/(n_samps * n_chans));
			fprintf(1, '%s: only read %d trials\n', data_file, n_trials_read);
		end
		is_mc_bin = true;
		if isfield(opt, 'reref_chan') && opt.reref_chan ~= 0
			dataX = bsxfun(@minus, data, data(:, opt.reref_chan));
			dataX(:, opt.reref_chan) = ...
				-data(:, opt.reref_chan);
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
			fprintf(1, '%s: error: no trial info provided\n', data_file);
			return
		end
		nr = size(trial_info, 1);
		if nr < n_trials % zero pad
			fprintf(1, ...
				'%s: error: trial info inconsistent with number of trials ', ...
				data_file);
			fprintf(1, 'zero padding trial info\n');
			nc = size(trial_info, 2);
			new_trial_info = zeros(n_trials, nc);
			new_trial_info(1:nr, :) = trial_info;
			trial_info = new_trial_info;
			clear new_trial_info
		elseif nr > n_trials
			fprintf(1, ...
				'%s: error: trial info inconsistent with number of trials ', ...
				data_file);
			fprintf(1, 'truncating trial info\n');
			trial_info = trial_info(1:n_trials, :);
		end
		trial_type = trial_info(:, opt.trial_col);
    
        
    % if working with a *.cnt.h1 file
	elseif strcmpi(data_file_type, 'hdf_binary')
        is_mc_bin=false;
        % first read in the HDF file
		h1_struct = read_hdf1_dataV7(data_file);
		if isempty(h1_struct) 
			fprintf(1, 'error reading %s\n', data_file);
			return;
		end
		X = h1_struct.data_struct.hdf1_cnt_data;
		if isempty(X)
			fprintf(1, 'no data in %s\n', data_file);
			return;
		else
			data = X';
			clear X
			% at this point data is n_total_samps X n_chans
		end
		% unload variables
		[trial_type, ev_nums] = check_hdf_trials(h1_struct);
        
	else
		fprintf(1, 'file type should be either "mc_binary" or "hdf_binary"\n');
		return;
    end
    
    % clip points off the end to mitigate filtering edge effects, always
    % comes out as the terminal 19 samples or so
    % (artifact from turning off amplifier most likely)
    data = clip_end_rawdata(data, 10000);

%%
    
    % if doing cleaning
    if opt.cleanset
        
        if isfield(opt,'hpfilt_bands')
            hpfilt_bands = opt.hpfilt_bands;
        else
            hpfilt_bands = [];
        end
        
        if isfield(opt, 'eog_rem') && ~isempty(opt.eog_rem) ...
                && isfield(opt, 'eog_chans') && ~isempty(opt.eog_chans)
            eog_data = data(:, opt.eog_chans);
            EOG=import_hdf_eeglab_inline(data_file, eog_data, h1_struct);
            EOG = clean_drifts(EOG, hpfilt_bands);
            eog_data = EOG.data';
            clear EOG
        else
            eog_data = [];
        end
        
        % cut out only the channels of interest
        if isfield(opt, 'chan_vec')
            dataC = data(:, opt.chan_vec);
            raw_size=size(dataC);
        else
            raw_size=size(data);
        end
        clear data
        
        %put the data into an EEGLAB dataset structure
        EEG=import_hdf_eeglab_inline(data_file, dataC, h1_struct);
        clear dataC
        
        if opt.vis
            %visualize with trimoutlier
            [EEG,h] = pop_trimOutlierPlot(EEG);
            print([subj_ckpath,name,'_1raw'],'-dpng');
            close(h)
        
            % high-pass filter at 1 Hz and visualize again (not used
            EEG2 = pop_eegfiltnew(EEG, [], 1, 1650, true, [], 0);
        
            %visualize with trimoutlier
            [EEG2,h] = pop_trimOutlierPlot(EEG2);
            print([subj_ckpath,name,'_2hpf'],'-dpng');
            close(h)
        end
        
        %save channel locations
        chan_locs=EEG.chanlocs;
        
        if isfield(opt,'burst')
            burst = opt.burst;
        else
            burst = 'off'; % off by default for event-related analyses
        end
        
        % filter and clean using ASR
        EEG = cleanset_eeglab(EEG, hpfilt_bands, burst);
        
        %interpolate any missing channels and record how many were bad
        n_interpchans=length(chan_locs) - size(EEG.data,1);
        EEG = pop_interp(EEG, chan_locs, 'spherical');
        
        if opt.vis
            %visualize with trimoutlier
            [EEG,h] = pop_trimOutlierPlot(EEG);
            print([subj_ckpath,name,'_3asr'],'-dpng');
            close(h)
        end
        
        % if desired, remove eyeblinks
        if isfield(opt, 'eog_rem') && ~isempty(opt.eog_rem)
            dataB = [EEG.data',eog_data];
            eog_chans = size(EEG.data,1)+1:size(EEG.data,1)+size(eog_data,2);
            clear EEG
            EEG=import_hdf_eeglab_inline(data_file, dataB, h1_struct);
            clear dataB
            if strcmpi(opt.eog_rem, 'lms')
                for mu = 6:8
                    try
                        EEG = pop_lms_regression( EEG, eog_chans, 3, 1.*10.^-mu, []);
                        fprintf(1, 'mu used was %d\n', mu);
                        break
                    catch
                        continue
                    end
                end
            elseif strcmpi(opt.eog_rem, 'crls')
                try
                    EEG = pop_crls_regression( EEG, eog_chans, 3, 0.9999, 0.01, []);
                catch
                    EEG = pop_scrls_regression( EEG, eog_chans, 3, 0.9999, 0.01, 50, []);
                end
            elseif strcmpi(opt.eog_rem, 'hinf')
                try
                    EEG = pop_hinfew_regression( EEG, eog_chans, 3, 0.005, 1e-05, 1.5, 0.99, []);
                catch e
                    fprintf(1, 'EOG failed: %s\n', e.message)
                end
            elseif strcmpi(opt.eog_rem, 'bss')
                try
                    EEG = pop_autobsseog( EEG, EEG.xmax, EEG.xmax, 'sobi', ...
                        {'eigratio', [1000000]}, 'eog_corr', {'range',[2  21]}, eog_chans);
                catch e
                    fprintf(1, 'EOG failed: %s\n', e.message)
                    if strcmpi(e.identifier, 'MATLAB:license:checkouterror')
                        return
                    end
                end
            end
            EEG.data(end-length(eog_chans)+1:end,:)=[];
            if strcmpi(opt.eog_rem, 'bss')
                EEG.chanlocs = chan_locs;
                start_thresh = 5;
                rej_elecs = zeros(13, 1);
                while length(rej_elecs) > 12
                    [EEG, rej_elecs] = pop_rejchan(EEG, 'elec',[1:61] ,'threshold', start_thresh, 'norm','on','measure','prob');
                    start_thresh = start_thresh + 1;
                end
                EEG = pop_interp(EEG, chan_locs, 'spherical');
                n_interpchans = n_interpchans + length(rej_elecs);
            end
        end
        
        dataF = EEG.data';
        clear EEG
        
        %check to make sure there are the same number of data points as before
        if any(size(dataF) ~= raw_size)
            fprintf(1, 'error: %s caused a problem with filtering and/or ASR\n', ...
                data_file);
            return
        end
        clear raw_size
    else
        
        % cut out only the channels of interest
        if isfield(opt, 'chan_vec')
            dataF = data(:, opt.chan_vec);
        else
            dataF = data;
        end
        clear data
    end
    
    n_chans = size(dataF, 2);
    
    % CSD transform the data if applicable
    if isfield(opt,'csd_G') && ~isempty(opt.csd_G) && isfield(opt,'csd_H') && ~isempty(opt.csd_H)
        %transpose in and out
        dataR = CSD(dataF', opt.csd_G, opt.csd_H, 1.0e-5, 10.0)'; %smooth lambda = -5, head size = 10 cm
        clear dataF
    else
        dataR = dataF;
        clear dataF
    end
    
    if opt.vis
        EEG=import_hdf_eeglab_inline(data_file,dataR,h1_struct);
        [EEG,h] = pop_trimOutlierPlot(EEG);
        print([subj_ckpath,name,'_4'],'-dpng');
        close(h)
        clear EEG
    end
    
    % memoize dataR here to save LOTS of time
    if isfield(opt, 'save_clean')
        if opt.save_clean
            if ~exist(cleanfile_dir, 'dir')
                mkdir(cleanfile_dir)
            end
            save(cleanfile_path, 'dataR', 'is_mc_bin', 'n_chans', 'trial_type', ...
                'n_interpchans', 'data_file', 'eog_data', 'ev_nums', '-v7.3')
        end
    end

else
    
    %load it and the associated h1
    load(cleanfile_path)   %should contain dataR, rate, is_mc_bin, n_chans, n_interpchans
    h1_struct = read_hdf1_dataV7(data_file);
end

% if not ASR-cleaned, filter data
if ~opt.cleanset
    
    % no channels will ever be interpolated if not ASR-cleaned
    n_interpchans = 0;
    
    % hp_cutoff-50 Hz using david's function
    hp_cutoff = mean(opt.hpfilt_bands);
    dataR = filter_data(dataR, hp_cutoff, 50, [], opt.rate);
    
    % eeglab way
    %EEG = import_hdf_eeglab_inline(data_file, dataR, h1_struct);
    %clear dataR
    %hp_cutoff = mean(opt.hpfilt_bands);
    %EEG = pop_eegfiltnew(EEG, hp_cutoff, []);
    %dataR = EEG.data';
end

% if not CSD-transformed, consider re-reference
if isfield(opt,'reref') && ~strcmpi(opt.reref,'none') && ...
    isfield(opt,'csd_G') && isempty(opt.csd_G) && isfield(opt,'csd_H') && isempty(opt.csd_H)
    %dataR = reference_data(dataR, opt.reref); %can be 'none' or 'average'
    
    % only average re-ref if requested
    if strcmpi(opt.reref, 'average')
        EEG = import_hdf_eeglab_inline(data_file, dataR, h1_struct);
        clear dataR
    
        EEG = pop_reref( EEG, []);
        dataR = EEG.data';
        clear EEG
    end
    
end

%%

% do wavelet calcuation
n_scales = numel(opt.wavelet_scales);
if isfield(opt, 'wavelet_cycles')
    if ~isempty(opt.wavelet_cycles)
        dataW = wavelet_calc2(dataR, opt.wavelet_scales, opt.wavelet_cycles);
    else
        dataW = wavelet_calc(dataR, opt.wavelet_scales);
    end
else
    dataW = wavelet_calc(dataR, opt.wavelet_scales);
end

if is_mc_bin
    dataW = reshape(dataW, [n_samps, n_trials, n_chans, n_scales]);
else

    %reshape wavelet data into trials
    [~, ~, dataW] = hdf_reshapeW(dataW, h1_struct, opt);
    %reshape cleaned and possibly CSD-transformed time-series into trials
    [~, n_samps, dataR, trial_type, ev_nums] = hdf_reshape(dataR, h1_struct, opt, trial_type, ev_nums);
    %reshape cleaned, pre-CSD time series into trials
    %[~, ~, dataF] = hdf_reshape(dataF, h1_struct, opt);

    if isempty(dataR)
        fprintf(1, 'invalid time data in %s\n', data_file);
        return
    end
    n_trials = size(dataR, 2);
    if n_trials < numel(trial_type)
        trial_type = trial_type(1:n_trials);
    end
end 

% dataR is n_samps X n_trials X n_chans

% determine the experiment for behavioral analysis (happens to be subdir1)
exp_name = subdir1;
exp_ind = find( strcmpi(exp_name, {opt.expstruct.name}) );

% if the experiment is recognizable
if ~isempty(exp_ind)

    % make an event table with first trial specified
    etable=h1_getbehav2(h1_struct, opt.expstruct(exp_ind));

    % if it was the gambling task, some custom behavior calculations
    % were coded
    if opt.trial_init_type==90 && length(opt.case_vec)==9
        [opt.trial_mat,opt.case_label,respRT] = ...
            behavinds_sog(etable);
    elseif opt.trial_init_type==90 && length(opt.case_vec)==8
        [opt.trial_mat,opt.case_label,respRT] = ...
            behavinds_sog2(etable);
    elseif opt.trial_init_type==90 && length(opt.case_vec)==4
        [~,~,respRT] = ...
            behavinds_sog(etable);
    %elseif opt.trial_init_type==90 && length(opt.case_vec)==2
    %    [~,~,respRT] = ...
    %        behavinds_sog(etable);
    end
else
% otherwise just make a generic event table
    etable=h1_getbehav(h1_struct);
end

n_cases = max(trial_type);
% if no custom trial indexings were done above
if ~isfield(opt, 'trial_mat')

    % make trial index based on type codes
    trial_mat = false(n_trials, n_cases);
    for m = 1:n_cases
        trial_mat(:, m) = trial_type == m;
    end
    if isfield(opt, 'case_vec') %in other words, if there is a
        % case vector specified, take only those, otherwise take all
        % cases
        trial_mat = trial_mat(:, opt.case_vec);
    end
    zero_trials = find(sum(trial_mat) == 0);
    if ~isempty(zero_trials)
        trial_mat(:, zero_trials) = [];
    end
else
    trial_mat = opt.trial_mat;
end
n_cases = size(trial_mat, 2);

% trial rejection

% first specify the channel / timepoints range to examine
if isfield(opt, 'artf_chan_vec')
    artf_chan_vec = opt.artf_chan_vec; %use specified channels if specified
else
    artf_chan_vec = 1:n_chans; %or all of them if not
end

if isfield(opt, 'artf_samp_vec')
    artf_samp_vec = opt.artf_samp_vec; %use certain samples if specified
elseif isfield(opt, 'artf_ms_vec')
    artf_samp_vec = round(opt.artf_ms_vec/1000*...
        opt.rate);
    artf_samp_vec=artf_samp_vec(1):artf_samp_vec(2); %or a ms range, better
end
if isfield(opt, 'prestim_ms')
    artf_samp_vec = artf_samp_vec + round(opt.prestim_ms/1000 * opt.rate) + 1;
elseif isfield(opt, 'epoch_lims')
    artf_samp_vec = artf_samp_vec + round(abs(opt.epoch_lims(1))/1000 * opt.rate) + 1;
else
    artf_samp_vec = 49:166;
end

% epoch rejection
trial_mat_orig=trial_mat;
if isfield(opt,'epochrej_method')
    
if strcmpi(opt.epochrej_method,'simple')
    trial_mat = check_artifact(dataR(artf_samp_vec, :, :), trial_mat_orig, ...
        artf_chan_vec, opt.thresh);
elseif strcmpi(opt.epochrej_method,'simple_bl')
    if isfield(opt, 'thresh_bl') && isfield(opt, 'artf_ms_vec')
        ms_vec = opt.artf_ms_vec(1):1000/opt.rate:opt.artf_ms_vec(2);
        ms_vec(end) = [];
        bl_pts = dsearchn(ms_vec', [-100, 0]');
        bl_vec = bl_pts(1):bl_pts(end);
    else
        bl_vec = 103:129;
    end
    trial_mat = check_artifact_baselined(dataR(artf_samp_vec, :, :), trial_mat_orig, ...
        artf_chan_vec, opt.thresh, bl_vec);
elseif strcmpi(opt.epochrej_method,'eeglab')
    for sdthresh=4 %:.25:6
    trial_mat = check_artifact_eeglab(data_file, dataR(artf_samp_vec, :, :), ...
        h1_struct, trial_mat_orig, sdthresh, opt.thresh, artf_chan_vec);
    prop_kept=sum(sum(trial_mat))/sum(sum(trial_mat_orig));
        if prop_kept >= .7
            fprintf(1, 'SD threshold used was %1.2f\n', sdthresh);
            break
        end
    end
elseif strcmpi(opt.epochrej_method,'both') && isfield(opt,'simple_thresh') ...
        && isfield(opt,'eeglab_thresh')
    if isfield(opt, 'csd_G') && ~isempty(opt.csd_G)
        uv_add_amt = 1;
    else
        uv_add_amt = 25;
    end
    trial_mat_both = repmat(trial_mat, [1 1 2]);
    prop_kept_simple = 0;
    while prop_kept_simple < 0.5
        trial_mat_both(:,:,1) = check_artifact(dataR(artf_samp_vec, :, :), trial_mat_orig, ...
            artf_chan_vec, opt.simple_thresh);
        prop_kept_simple = sum(sum(trial_mat_both(:,:,1)))/sum(sum(trial_mat_orig));
        opt.simple_thresh = opt.simple_thresh + uv_add_amt;
    end
    fprintf(1, 'Simple uV threshold used was %1.2f\n', opt.simple_thresh - uv_add_amt);
    prop_kept_eeglab = 0;
    start_sd_thresh = 4;
    while prop_kept_eeglab < 0.5
        trial_mat_both(:,:,2) = check_artifact_eeglab(data_file, dataR(artf_samp_vec, :, :), ...
            h1_struct, trial_mat_orig, start_sd_thresh, opt.eeglab_thresh, artf_chan_vec);
        prop_kept_eeglab = sum(sum(trial_mat_both(:,:,2)))/sum(sum(trial_mat_orig));
        opt.eeglab_thresh = opt.eeglab_thresh + uv_add_amt;
        start_sd_thresh = start_sd_thresh + 0.25;
    end
    fprintf(1, 'EEGLAB uV threshold used was %1.2f\n', opt.eeglab_thresh - uv_add_amt);
    trial_mat = squeeze(all(trial_mat_both, 3));
    prop_kept=sum(sum(trial_mat))/sum(sum(trial_mat_orig));
    if prop_kept<=.4
        trial_mat = squeeze(any(trial_mat_both, 3)); %experimental
    end
end

else
    fprintf(1, 'trial rejection method incorrectly or not specified for %s\n', ...
        data_file);
    return
end

% downsample TF results
ds_rate=opt.tf_timedownsamp_ratio;
n_tfsamps = ceil( n_samps / ds_rate );

%put the data into nice matrices
mean_data = zeros(n_samps, n_chans, n_cases);
% wavelet_evk = zeros(n_tfsamps, n_chans, n_scales, n_cases);
wavelet_evknorm = zeros(n_tfsamps, n_chans, n_scales, n_cases);
% wavelet_tot = zeros(n_tfsamps, n_chans, n_scales, n_cases);
wavelet_totpow = zeros(n_tfsamps, n_chans, n_scales, n_cases);

for m = 1:n_cases
    mean_data(:, :, m) = squeeze(mean(dataR(:, trial_mat(:, m), :), 2)); %erps
    %wavelet_evk(:, :, :, m) = ...
    %    squeeze(mean(dataW(1:ds_rate:end, trial_mat(:, m), :, :), 2)); %weighted phase
    wavelet_evknorm(:,:,:,m) = squeeze(mean( ... %normalized phase
        dataW(1:ds_rate:end, trial_mat(:, m), :, :) ./ ...
        abs(dataW(1:ds_rate:end, trial_mat(:, m), :, :)), 2));
    %wavelet_tot(:, :, :, m) = ... %amplitude
    %    squeeze(mean(abs(dataW(1:ds_rate:end, trial_mat(:, m), :, :)), 2));
    %wavelet_totpow(:, :, :, m) = ... %power
    %    squeeze(mean(abs(dataW(1:ds_rate:end, trial_mat(:, m), :, :)).^2, 2)); %square to give power (uV^2)
    wavelet_totpow(:, :, :, m) = ... %power
        squeeze(mean(dataW(1:ds_rate:end, trial_mat(:, m), :, :).* ...
        conj(dataW(1:ds_rate:end, trial_mat(:, m), :, :)), 2)); % this is a faster abs().^2
end

%if any(ismember(opt.measures,'wave_evk')) %deprecated
%    Y.wave_evk = wavelet_evk;
%end
if any(ismember(opt.measures,'wave_evknorm'))
    Y.wave_evknorm = single(wavelet_evknorm);
end
%if any(ismember(opt.measures,'wave_tot'))
%    Y.wave_tot = wavelet_tot;
%end
if any(ismember(opt.measures,'wave_totpow'))
    Y.wave_totpow = single(wavelet_totpow);
end
if any(ismember(opt.measures,'erptrial')) && exist('dataR','var')
    %[~, ~, dataF] = hdf_reshape(dataF, h1_struct, opt);
    Y.erptrial = single(dataR(:, any(trial_mat,2), :)); %accept only good trials
    Y.erptrial_log = trial_mat(any(trial_mat,2),:);
    Y.erptrial_evnums = ev_nums(any(trial_mat,2));
end

Y.erp = mean_data;
Y.trials = trial_mat;
Y.trials_both = trial_mat_both;
Y.n_trials = sum(trial_mat);
Y.n_interpchans = n_interpchans;
Y.opt = opt; % save opt now -- good idea!


% Y.etable = etable;
% interpret behavioral data
Y.beh = beh_etable( etable, opt.expstruct(exp_ind) );

%behavioral data saving
if isfield(opt,'getbehav')
    if opt.getbehav && opt.trial_init_type==90 && n_cases==4

    [avgbet_prevoutcome,stdbet_prevoutcome,crit,crit_prevoutcome,avgbet,resp_mat, ...
        prevoutcome_mat,avgbet2_prevoutcome,stdbet2_prevoutcome] = behav_sog(etable);

    behav_data=v2struct(avgbet_prevoutcome,stdbet_prevoutcome,crit,crit_prevoutcome, ...
        avgbet,resp_mat, prevoutcome_mat,avgbet2_prevoutcome,stdbet2_prevoutcome,respRT);

    Y.behav = behav_data;

    end
end
%Y.wavetrial = dataW;

if isfield(opt, 'n_cohperms') && isfield(opt, 'rng_seed')
    n_cohperms = opt.n_cohperms;
    rng_seed = opt.rng_seed;
else
    n_cohperms = 0;
    rng_seed = 0;
end

% calculate channel coherence
if isfield(opt,'coherence_pairs')
    coherence_pairs = opt.coherence_pairs;
    if isfield(opt, 'coherence_type')
        coherence_type = opt.coherence_type;
    else
        coherence_type = 'weighted';
    end
    if isfield(opt, 'chan_vec') % check for consistency with pairs
        [U, ~, Un] = unique(opt.coherence_pairs);
        if isequal(U, opt.chan_vec)
            coherence_pairs = reshape(Un, size(opt.coherence_pairs));
        elseif all(ismember(U, 1:n_chans))
            coherence_pairs = opt.coherence_pairs;
        else
            fprintf(1, ...
                'Inconsistent specification of chan_vec and coherence_pairs\n');
            coherence_pairs = [];
        end
    end
    if ~isempty(coherence_pairs)
        n_coh_pairs = size(opt.coherence_pairs, 1);
        if n_cohperms > 0
            coh_results{1} = zeros(n_tfsamps, n_scales, n_cases, n_coh_pairs);
            coh_results{2} = zeros(n_tfsamps, n_scales, n_cases, n_coh_pairs);
            for m = 1:n_coh_pairs
                [coh_results{1}(:, :, :, m),  coh_results{2}(:, :, :, m)] = ...
                    pairwise_coherence_calc(dataW(1:ds_rate:end, :, coherence_pairs(m, :), :), ...
                        trial_mat, n_cohperms, rng_seed, coherence_type);
            end
        else
            coh_results = zeros(n_tfsamps, n_scales, n_cases, n_coh_pairs);
            for m = 1:n_coh_pairs
                coh_results(:, :, :, m) = ...
                    pairwise_coherence_calc(dataW(1:ds_rate:end, :, coherence_pairs(m, :), :), ...
                        trial_mat, n_cohperms, rng_seed, coherence_type);
            end		
        end
    else
        coh_results = [];
    end
    if any(ismember(opt.measures,'coh'))
        Y.coh = single(coh_results);
    end
end

%end of function
end
