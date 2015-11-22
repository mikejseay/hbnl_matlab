function [Y,name] = beh_batch(data_file,data_file_type,param_struct)

%get the subject id
[~,name,~]=fileparts(data_file);

% stage 1
	Y = [];
	% get parameter from param_struct
    
    %check for existence of cleaned raw H1 data as "dataR" and any other needed vars
% build a suffix based on presence of cleaning and CSD
check_suffix='';
if isfield(param_struct,'cleanset')
    if param_struct.cleanset
        check_suffix=[check_suffix,'_clean'];
    end
end
if isfield(param_struct,'csd_G')
    check_suffix=[check_suffix,'_CSD'];
end
    % if the cleaned data corresponding to the session exists already
if isfield(param_struct,'cleandir')    
    cleanfile_path=[param_struct.cleandir,'/',name,check_suffix,'.mat'];
else
    cleanfile_path=[];
    fprintf('no directory for cleaned h1 recordings was specified');
end
% if the data has already been cleaned
if exist(cleanfile_path,'file')
    %load it
    load(cleanfile_path)   %should contain dataR and h1_struct
else % else calculate it
        
	if exist(data_file, 'file') ~= 2
		fprintf(2, 'error: %s does not exist\n', data_file);
		return
    end	    
    
	if strcmpi(data_file_type, 'mc_binary')
        is_mc_bin=true;
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
        is_mc_bin=false;
		h1_struct = read_hdf1_dataV7(data_file);
		if isempty(h1_struct) 
			fprintf(2, 'error reading %s\n', data_file);
			return;
		end
		
		% unload variables
		trial_type = check_hdf_trials(h1_struct);
	else
		fprintf(2, 'only mc_binary files can be used at this point\n');
		return;
	end
end
        
    if isfield(param_struct, 'trial_init_type')
        etable=h1_getbehav(h1_struct,param_struct.trial_init_type);
        if param_struct.trial_init_type==90 && length(param_struct.case_vec)==9
            [param_struct.trial_mat,param_struct.case_label,respRT] = ...
                behavinds_sog(etable);
        elseif param_struct.trial_init_type==90 && length(param_struct.case_vec)==8
            [param_struct.trial_mat,param_struct.case_label,respRT] = ...
                behavinds_sog2(etable);
        elseif param_struct.trial_init_type==90
            [~,~,respRT] = ...
                behavinds_sog(etable);
        end
    else
        etable=h1_getbehav(h1_struct);
    end
	
	%check to see if last trial was removed and reflect that change in the
	%trial_mat
    %if length(trial_indexW) ~= size(trial_mat,1)
    %    trial_mat = trial_mat(1:(end - 1), :);
    %end
    
% end of stage 1	
% filter and calculate CWT

	% eliminate first and last trials
	
    Y.etable = etable;
    %behavioral data saving
    if isfield(param_struct,'getbehav')
        if param_struct.getbehav && param_struct.trial_init_type==90
            
        [avgbet_prevoutcome,stdbet_prevoutcome,crit,crit_prevoutcome,avgbet,resp_mat, ...
 prevoutcome_mat,avgbet2_prevoutcome,stdbet2_prevoutcome] = behav_sog(etable);

        behav_data=v2struct(avgbet_prevoutcome,stdbet_prevoutcome,crit,crit_prevoutcome, ...
            avgbet,resp_mat, prevoutcome_mat,avgbet2_prevoutcome,stdbet2_prevoutcome,respRT);
        
        Y.behav_data = behav_data;
        
        end
    end

