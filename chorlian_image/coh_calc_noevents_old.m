function Y = coh_calc_noevents(data_file, opt)
% pre-processes EEG data, applies wavelet transform, and calculates
% coherence (ITC and ERPCOH).

% get the subject id
[~,name,~]=fileparts(data_file);

% initialize output structure
Y = [];
    
% check for existence of needed data (cleaned and/or CSD-transformed)
% build a suffix based on presence of cleaning and CSD
subdir1 = name(1:3);
check_suffix='';
if isfield(opt,'cleanset')
    if opt.cleanset
        check_suffix=[check_suffix,'_clean'];
        subdir2 = 'clean';
    else
        check_suffix=[check_suffix,'_noclean'];
        subdir2 = 'noclean';
    end
end
if isfield(opt,'csd_G') && ~isempty(opt.csd_G) && isfield(opt,'csd_H') && ~isempty(opt.csd_H)
    check_suffix=[check_suffix,'_CSD'];
    subdir2 = 'clean_CSD';
end
clean_name = [name, check_suffix, '.mat'];
if isfield(opt,'cleandir')    
    cleanfile_dir = strjoin( {opt.cleandir, subdir1, subdir2}, filesep);
    cleanfile_path=fullfile(opt.cleandir, subdir1, subdir2, clean_name);
else
    cleanfile_dir = [];
    cleanfile_path = [];
    fprintf('no directory for cleaned h1 recordings was specified\n');
end

% if the data has already been cleaned
if ~exist(cleanfile_path,'file')
%%    
    
    % if the file doesn't exist, error
    if exist(data_file, 'file') ~= 2
        fprintf(2, 'error: %s does not exist\n', data_file);
        return
    end
    
    is_mc_bin=false;

    EEG = pop_loadcnt(data_file, 'dataformat', 'auto');
    %rate = EEG.srate;
    
    EEG = pop_resample(EEG, opt.rate);
    data = EEG.data';

    %n_samps = round(range(opt.epoch_lims)*rate/1000);
    %[wavelet_scales, cycles] = findCWTscales(n_samps,rate,opt.freq_lims,opt.padratio);
    
    clear EEG
        	
    % cut out only the channels of interest
	if isfield(opt, 'chan_vec')
        %try 
            data = data(:, opt.chan_vec);
        %catch
        %    warning('File lacks channels specified');
        %end
	end
	n_chans = size(data, 2);
    
    % clip points off the end to mitigate filtering edge effects, always
    % comes out as the terminal 19 samples or so
    % (artifact from turning off amplifier most likely)
    % data = clip_end_rawdata(data,10000);
    
    raw_size=size(data);
    
    % if doing cleaning using ASR
    if opt.cleanset
        
        %put the data into an EEGLAB dataset structure
        EEG=import_rawcnt_eeglab_inline(data, opt.rate);
        clear data
        
        %save channel locations
        chan_locs=EEG.chanlocs;
        
        if isfield(opt,'hpfilt_bands')
            hpfilt_bands = opt.hpfilt_bands;
        else
            hpfilt_bands = [];
        end
        
        if isfield(opt,'burst')
            burst = opt.burst;
        else
            burst = 'off'; % for consistency's sake
        end
        
        % filter and clean using ASR
        EEG = cleanset_eeglab(EEG, hpfilt_bands, burst);
        
        %interpolate any missing channels and record how many were bad
        n_interpchans=length(chan_locs) - size(EEG.data,1);
        EEG = pop_interp(EEG, chan_locs, 'spherical');
        
        dataF = EEG.data';
        clear EEG
    else
        dataF = data;
        clear data
    end
    
    %check to make sure there are the same number of data points as before
    if any(size(dataF) ~= raw_size)
        fprintf(2, 'error: %s caused a problem with filtering and/or ASR\n', ...
            data_file);
        return
    end
    clear raw_size
    
    % CSD transform the data if applicable
    if isfield(opt,'csd_G') && ~isempty(opt.csd_G) && isfield(opt,'csd_H') && ~isempty(opt.csd_H)
        %transpose in and out
        dataR = CSD(dataF', opt.csd_G, opt.csd_H, 1.0e-5, 10.0)'; %smooth lambda = -5, head size = 10 cm
    else
        dataR = dataF;
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

        if strcmpi(opt.reref, 'average')
            EEG = import_hdf_eeglab_inline(data_file, dataR, opt.rate);
            clear dataR

            EEG = pop_reref( EEG, []);
            dataR = EEG.data';
        end

    end
    
    % memoize dataR here to save LOTS of time
    if isfield(opt, 'save_clean')
        if opt.save_clean
            if ~exist(cleanfile_dir, 'dir')
                mkdir(cleanfile_dir)
            end
            save(cleanfile_path, 'dataR','is_mc_bin','n_chans','n_interpchans')
        end
    end

else
    
    %load it
    load(cleanfile_path)   %should contain dataR, is_mc_bin, n_chans, n_interpchans

end

%%

% do wavelet calcuation
n_scales = length(opt.wavelet_scales);
if isfield(opt, 'wavelet_cycles')
    if ~isempty(opt.wavelet_cycles)
        dataW = wavelet_calc2(dataR, opt.wavelet_scales, opt.wavelet_cycles);
    else
        dataW = wavelet_calc(dataR, opt.wavelet_scales);
    end
else
    dataW = wavelet_calc(dataR, opt.wavelet_scales);
end

%reshape into n_samps long segments (overlapping most likely)
%n_samps = opt.n_samps;
[~, dataW] = noevents_reshapeW(dataW, opt.n_samps);
[~, dataR] = noevents_reshape(dataR, opt.n_samps);

%get the number of segments made
n_trials = size(dataR, 2);

% dataR is n_samps X n_trials X n_chans
    
trial_mat = true(n_trials,1);
n_cases = 1;
    
% trial rejection

% first specify the channel / timepoints range to examine
if isfield(opt, 'artf_chan_vec')
    artf_chan_vec = opt.artf_chan_vec; %use specified channels if specified
else
    artf_chan_vec = 1:n_chans; %or all of them if not
end

%use all samples by default
artf_samp_vec=1:opt.n_samps;
    
% epoch rejection
trial_mat_orig=trial_mat;
h1_struct = struct();
h1_struct.experiment_struct.rate = opt.rate; %hack to make eeglab rejection work below

if isfield(opt,'epochrej_method')
    
if strcmpi(opt.epochrej_method,'simple')
    trial_mat = check_artifact(dataR(artf_samp_vec, :, :), trial_mat_orig, ...
        artf_chan_vec, opt.thresh);           
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
    trial_mat_both = repmat(trial_mat, [1 1 2]);
    trial_mat_both(:,:,1) = check_artifact(dataR(artf_samp_vec, :, :), trial_mat_orig, ...
        artf_chan_vec, opt.simple_thresh);
    trial_mat_both(:,:,2) = check_artifact_eeglab(data_file, dataR(artf_samp_vec, :, :), ...
        h1_struct, trial_mat_orig, 4.5, opt.eeglab_thresh, artf_chan_vec);
    trial_mat = squeeze(all(trial_mat_both, 3));
    prop_kept=sum(sum(trial_mat))/sum(sum(trial_mat_orig));
    if prop_kept<=.7
        trial_mat = squeeze(any(trial_mat_both, 3)); %experimental
    end
end

else
    fprintf(1, 'trial rejection method incorrectly or not specified for %s\n', ...
        data_file);
    return
end

% downsample TF results
%ds_rate=opt.tf_timedownsamp_ratio;
%n_tfsamps = ceil( n_samps / ds_rate );

%put the data into nice matrices
%wavelet_tot = zeros(n_chans, n_scales, n_cases);
wavelet_totpow = zeros(n_chans, n_scales, n_cases);

for m = 1:n_cases
    %wavelet_tot(:, :, m) = ... %amplitude
    %    squeeze(mean(mean(abs(dataW(:, trial_mat(:, m), :, :)), 1), 2));
    wavelet_totpow(:, :, m) = ... %power
        squeeze(mean(mean(dataW(:, trial_mat(:, m), :, :).* ...
        conj(dataW(:, trial_mat(:, m), :, :)),1 ), 2)); % this is a faster abs().^2
end

%if any(ismember(opt.measures,'wave_tot'))
%    Y.wave_tot = wavelet_tot;
%end
if any(ismember(opt.measures,'wave_totpow'))
    Y.wave_totpow = wavelet_totpow; % output is n_chans x n_scales
end

Y.trials = trial_mat;
Y.n_trials = sum(trial_mat);
Y.n_interpchans = n_interpchans;

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
            fprintf(2, ...
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
            coh_results = zeros(n_scales, n_cases, n_coh_pairs);
            for m = 1:n_coh_pairs
                coh_results(:, :, m) = ...
                    pairwise_coherence_calc_noevents(dataW(:, :, coherence_pairs(m, :), :), ...
                        trial_mat, n_cohperms, rng_seed, coherence_type);
            end		
        end
    else
        coh_results = [];
    end
    if any(ismember(opt.measures,'coh'))
        Y.coh = squeeze(coh_results)'; % output is n_pairs x n_scales
    end
end

%end of function
end
