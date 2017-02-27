function calc_hdf1_st_v40_ms(cnthdf1_filename, working_directory, the_case, chan_list_type,baseline_type,st_baseline_time_min_ms,st_baseline_time_max_ms, min_trials, max_trials, new_pre_stim_time_ms, new_threshold_value,new_threshold_min_time_ms, new_threshold_max_time_ms, new_hi_pass, new_lo_pass, new_threshold_electrodes) 
%
% function calc_hdf1_st_v4.0(cnthdf1_filename, [working_directory], [the_case],
% [chan_list_type], [baseline_type], [st_baseline_time_min_ms], 
% [st_baseline_time_max_ms], [min_trials], [max_trials], [new_pre_stim_time_ms], 
% [new_threshold_value], [new_threshold_min_time_ms], [new_threshold_max_time_ms], [new_hi_pass], [new_lo_pass],
% [new_threshold_electrodes]) 
%
%
% adding of new_hi_pass and new_lo_pass variables, Niklas (2008-07-21)
% adding of new_threshold_electrodes variable, Niklas (2008-07-21)
% removed the reading of /export/age_gender_file/exp_session_subject_data.csv into h1-file (2010-11-23)

 if (nargin < 16)
    new_threshold_electrodes = [];
 end
 if (nargin < 15)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
 end
 if (nargin < 14)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
 end
 if (nargin < 13)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
 end 
 if (nargin < 12)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
    new_threshold_min_time_ms = 0;
 end 
 if (nargin < 11)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
    new_threshold_min_time_ms = 0;
    new_threshold_value = 0;
 end 
 if (nargin < 10)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
    new_threshold_min_time_ms = 0;
    new_threshold_value = 0;
    new_pre_stim_time_ms = 0;  
 end 
 if (nargin < 9)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
    new_threshold_min_time_ms = 0;
    new_threshold_value = 0;
    new_pre_stim_time_ms = 0;  
    max_trials = 35;
 end 
 if (nargin < 8)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
    new_threshold_min_time_ms = 0;
    new_threshold_value = 0;
    new_pre_stim_time_ms = 0;  
    max_trials = 35;
    min_trials = 1;
 end 
 if (nargin < 7)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
    new_threshold_min_time_ms = 0;
    new_threshold_value = 0;
    new_pre_stim_time_ms = 0;  
    max_trials = 35;
    min_trials = 1;
    st_baseline_time_max_ms = 0;
 end 
 if (nargin < 6)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
    new_threshold_min_time_ms = 0;
    new_threshold_value = 0;
    new_pre_stim_time_ms = 0;  
    max_trials = 35;
    min_trials = 1;
    st_baseline_time_max_ms = 0;
    st_baseline_time_min_ms = 0;
 end 
 if (nargin < 5)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
    new_threshold_min_time_ms = 0;
    new_threshold_value = 0;
    new_pre_stim_time_ms = 0;  
    max_trials = 35;
    min_trials = 1;
    st_baseline_time_max_ms = 0;
    st_baseline_time_min_ms = 0;
    baseline_type = 1;
 end 
 if (nargin < 4)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
    new_threshold_min_time_ms = 0;
    new_threshold_value = 0;
    new_pre_stim_time_ms = 0;  
    max_trials = 35;
    min_trials = 1;
    st_baseline_time_max_ms = 0;
    st_baseline_time_min_ms = 0;
    baseline_type = 1;
    chan_list_type = 1;
 end 
 if (nargin < 3)
    new_threshold_electrodes = [];
    new_lo_pass = 55.0;
    new_hi_pass = 0.05;
    new_threshold_max_time_ms = 0;
    new_threshold_min_time_ms = 0;
    new_threshold_value = 0;
    new_pre_stim_time_ms = 0;  
    max_trials = 35;
    min_trials = 1;
    st_baseline_time_max_ms = 0;
    st_baseline_time_min_ms = 0;
    baseline_type = 1;
    chan_list_type = 1;
    the_case = 1;
 end 
 if (nargin < 2)
   new_threshold_electrodes = [];
   new_lo_pass = 55.0;
   new_hi_pass = 0.05;
   new_threshold_max_time_ms = 0;
   new_threshold_min_time_ms = 0;
   new_threshold_value = 0;
   new_pre_stim_time_ms = 0;  
   max_trials = 35;
   min_trials = 1;
   st_baseline_time_max_ms = 0;
   st_baseline_time_min_ms = 0;
   baseline_type = 1;
   chan_list_type = 1;
   the_case = 1;
   
end 
if (new_pre_stim_time_ms < 0) %in case someone used the minus sign with this value
    new_pre_stim_time_ms = abs(new_pre_stim_time_ms);
end

if (isunix == true)
    working_directory = [working_directory, '/'];
else
    working_directory = [working_directory, '\'];
end  
 
slash = findstr(cnthdf1_filename, '/');
if (length(slash) == 0)
    slash = 0;
end

sub_id     = cnthdf1_filename(slash(length(slash))+1:length(cnthdf1_filename));
underscore = findstr(sub_id, '_');	
subject_id = sub_id(underscore(3)+1:underscore(4)-1);
subject_session_run = sub_id(underscore(2)+1:underscore(2)+2);

% changed program from David to run Kevin's old matlab programs with matlab2010.
h1_struct = read_hdf1_dataV7(cnthdf1_filename);

% change post stim time to 1125 (always) - mike, 12/5/16
h1_struct.experiment_struct.post_stim_time_ms = 1125;

if (h1_struct.file_struct.data_type ~= 0)  
    return
end
if (new_threshold_value == -1) 
    h1_struct.experiment_struct.apply_thresholding = 0;
end    
if (new_threshold_value > 0) 
    h1_struct.experiment_struct.apply_thresholding = 1; 
    h1_struct.experiment_struct.threshold_value = new_threshold_value;
end

if (new_threshold_min_time_ms < 0) 
    h1_struct.experiment_struct.threshold_time_min_ms = new_threshold_min_time_ms; 
end

if (new_threshold_max_time_ms > 0) 
    h1_struct.experiment_struct.threshold_time_max_ms = new_threshold_max_time_ms; 
end

%if (st_baseline_time_max_ms == 0 & st_baseline_time_min_ms == 0)
%    st_baseline_time_min_ms = h1_struct.experiment_struct.baseline_time_min_ms;
%    st_baseline_time_max_ms = h1_struct.experiment_struct.baseline_time_max_ms;
%else
%    h1_struct.experiment_struct.baseline_time_min_ms = st_baseline_time_min_ms;
%    h1_struct.experiment_struct.baseline_time_max_ms = st_baseline_time_max_ms;
%end  
  
if (st_baseline_time_min_ms == 0)
    st_baseline_time_min_ms = h1_struct.experiment_struct.baseline_time_min_ms;
elseif ( st_baseline_time_min_ms < -h1_struct.experiment_struct.pre_stim_time_ms )
    st_baseline_time_min_ms = -h1_struct.experiment_struct.pre_stim_time_ms;
else
    h1_struct.experiment_struct.baseline_time_min_ms = st_baseline_time_min_ms;   
end    

if (st_baseline_time_max_ms == 0)
    st_baseline_time_max_ms = h1_struct.experiment_struct.baseline_time_max_ms;
else
    h1_struct.experiment_struct.baseline_time_max_ms = st_baseline_time_max_ms;   
end    

if (new_pre_stim_time_ms ~= 0)
    h1_struct.experiment_struct.pre_stim_time_ms = new_pre_stim_time_ms;
end  
pre_stim_time_ms = h1_struct.experiment_struct.pre_stim_time_ms;


%%%%% opening txt-file to export data of the mat-file creation process   
%txt_file_name = ['calc_hdf1_data_', subject_id,'_',h1_struct.experiment_struct.exp_name,'_',subject_session_run,'.log'];
txt_file_dir = '/active_projects/ERO_scripts/matlab_logs_cnth1st/';
txt_file_fname = [sub_id,'.st.mat.v4.0.log'];
txt_file_name = fullfile(txt_file_dir, txt_file_fname);
txt_exp       = cnthdf1_filename(slash(length(slash))+1:length(cnthdf1_filename)-25);

fid = fopen(txt_file_name, 'wt');
	fprintf(fid, 'Information about the mat-file creation process: \n' );
	fprintf(fid, '\n');
	fprintf(fid, 'Subject ID:           %s \n', subject_id);
    fprintf(fid, 'Experiment name:      %s \n', txt_exp);
    fprintf(fid, 'Experiment session:   %s \n', subject_session_run);
    %fprintf(fid, 'Experiment condition: %s \n', num2str(the_case));
	fprintf(fid, '\n');
   	fprintf(fid, 'Pre stimulus time :  %6.1f ms \n', h1_struct.experiment_struct.pre_stim_time_ms);
    fprintf(fid, 'Post stimulus time:  %6.1f ms \n', h1_struct.experiment_struct.post_stim_time_ms);
	fprintf(fid, '\n');
   	fprintf(fid, 'Threshold min time: %6.1f ms \n', h1_struct.experiment_struct.threshold_time_min_ms);
   	fprintf(fid, 'Threshold max time: %6.1f ms \n', h1_struct.experiment_struct.threshold_time_max_ms);
   	fprintf(fid, 'Threshold value   : %4.0f uV \n', h1_struct.experiment_struct.threshold_value);
	fprintf(fid, '\n');
    fprintf(fid, 'Lowpass filter before error checking:  %5.2f Hz \n', new_lo_pass);
   	fprintf(fid, 'Highpass filter before error checking: %5.2f Hz \n', new_hi_pass);
	fprintf(fid, '\n');
    fprintf(fid, 'Baseline min time: %6.1f ms \n', h1_struct.experiment_struct.baseline_time_min_ms);
    fprintf(fid, 'Baseline max time: %6.1f ms \n', h1_struct.experiment_struct.baseline_time_max_ms);
    fprintf(fid, '\n');

    h1_struct = create_cnthdf1_baseline_epoch_data_v4_ms(h1_struct, txt_file_name, new_hi_pass, new_lo_pass, new_threshold_electrodes);
    h1_struct = determine_correct_and_accepted_trials_v40(h1_struct, txt_file_name);

%fclose(fid);

h1_struct = create_epoch_average_data(h1_struct);

clear h1_struct.data_struct.hdf1_ind_st_data;

n_trials  = size(h1_struct.data_struct.hdf1_epoch_data,1);
n_samples = size(h1_struct.data_struct.hdf1_epoch_data,3);
n_cases   = h1_struct.experiment_struct.n_cases;
rate      = h1_struct.experiment_struct.rate;

case_idx      = get_case_idx(h1_struct.case_struct.case_num, the_case);
exp_case_type = lower(strtrim(h1_struct.case_struct.case_type{case_idx}));

%fid = fopen(txt_file_name, 'a');
    	fprintf(fid, '\n');
   	fprintf(fid, 'Case numbers               :');
	for i = 1:n_cases
   	fprintf(fid, '%3.0f ', h1_struct.case_struct.case_num(i));
	end
    	fprintf(fid, '\n');
   	fprintf(fid, 'Case names                 :');
	for i = 1:n_cases
   		if (length(h1_struct.case_struct.case_type{i}) == 1) 
			fprintf(fid, '  %s ', h1_struct.case_struct.case_type{i});
		elseif (length(h1_struct.case_struct.case_type{i}) == 2)
			fprintf(fid, ' %s ',  h1_struct.case_struct.case_type{i});
		else
			fprintf(fid, '%s ',   h1_struct.case_struct.case_type{i});
		end
	end
    	fprintf(fid, '\n');
   	fprintf(fid, 'Number of trials (all)     :');
	for i = 1:n_cases
   		fprintf(fid, '%3.0f ', h1_struct.case_struct.n_trials(i));
	end
    	fprintf(fid, '\n');
   	fprintf(fid, 'Number of trials (response):');
	for i = 1:n_cases
   		fprintf(fid, '%3.0f ', h1_struct.case_struct.n_responses(i));
	end
    	fprintf(fid, '\n');
   	fprintf(fid, 'Number of trials (omitted) :');
	for i = 1:n_cases
   		fprintf(fid, '%3.0f ', h1_struct.case_struct.n_resp_omitted(i));
	end
    	fprintf(fid, '\n');
   	fprintf(fid, 'Number of trials (errors)  :');
	for i = 1:n_cases
   		fprintf(fid, '%3.0f ', h1_struct.case_struct.n_resp_errors(i));
	end
    	fprintf(fid, '\n');
   	fprintf(fid, 'Number of trials (accepted):');
	for i = 1:n_cases
   		fprintf(fid, '%3.0f ', h1_struct.case_struct.n_trials_accepted(i));
	end
    	fprintf(fid, '\n');
	fprintf(fid, '\n');

fclose(fid);

if (case_idx == 0)
   return 
end    

case_accepted_idx = get_accepted_case_trial_idxs(h1_struct, the_case);
case_accepted_idx = case_accepted_idx(randperm(length(case_accepted_idx)));

if (length(case_accepted_idx) == 0)
    return
end

if (length(case_accepted_idx) < min_trials)
    return
end


% if (h1_struct.experiment_struct.exp_name == 'cpt')
 if ( txt_exp == 'cpt')
     if (exp_case_type == 't' )
         exp_case_type = 'cg';
     end
     if (exp_case_type == 'ng' )
        exp_case_type = 'un';
     end
     if (exp_case_type == 'n' )
         exp_case_type = 'dn';
     end
 end
% if (h1_struct.experiment_struct.exp_name == 'vp3')
 if (txt_exp == 'vp3')
     if (exp_case_type == 't' )
         exp_case_type = 'tt';
     end
     if (exp_case_type == 'n' )
         exp_case_type = 'nv';
     end
 end

fname = [working_directory, sub_id, '.', exp_case_type, '.st.mat'];

if (exist(fname, 'file') ~= 0) % st file already exits
   return 
end    

num_trials = min([length(case_accepted_idx) max_trials]);

[channel_list, n_chans] = get_channel_list(chan_list_type);
channel_indices = get_channel_indices(h1_struct.run_struct.channel_label, channel_list);

%if (new_hi_pass > 0.0)
   Wn_hi = new_hi_pass / (rate/2);
   [bhi,ahi] = cheby1(3,0.5,Wn_hi,'high');
%end    

for i=1:n_chans
    channel_idx = channel_indices(i);
    channel_name = get_channel_name(h1_struct.run_struct.channel_label, channel_idx);
    
    sprintf('%s', [num2str(i), ': ', channel_name{:}])
    
    if (channel_idx == 0)
       continue 
    end    
	
    h1_struct.data_struct.hdf1_avg_st_data(i,1:size(h1_struct.data_struct.hdf1_avg_data,3)) = squeeze(h1_struct.data_struct.hdf1_avg_data(case_idx,channel_idx,:));    

    test_st = false;
    if (test_st == true)
        case_accepted_idx = zeros(1,size(case_accepted_idx,2)) + case_accepted_idx(1);
        h1_struct.data_struct.hdf1_avg_st_data(i,1:size(h1_struct.data_struct.hdf1_epoch_data,3)) = squeeze(h1_struct.data_struct.hdf1_epoch_data(case_accepted_idx(1),channel_idx,:));    
    end
    
    trial_data = squeeze(h1_struct.data_struct.hdf1_epoch_data(case_accepted_idx,channel_idx,:));
     
	mean_data = mean(trial_data, 1);
	%if (new_hi_pass > 0)
		mean_data = filtfilt(bhi,ahi,mean_data); 
	%end  
	[Savg,t,f] = stockwell_transform(mean_data, 0, length(mean_data)/2, 1/rate); 
	%abs_Savg   = abs(Savg);
    %angle_Savg = angle(Savg);
    
    S_trial_vec = zeros(num_trials, size(Savg, 1), size(Savg, 2));
    for j=1:num_trials
       trial_vector = trial_data(j,:);
		%if (new_hi_pass > 0)
			trial_vector = filtfilt(bhi,ahi,trial_vector); 
		%end
       [S,t,f] = stockwell_transform(trial_vector, 0, length(trial_vector)/2, 1/rate); 
       
       S_trial_vec(j, :, :) = S;
       %abs_S   = abs(S);
       %angle_S = angle(S);
       
       % base_S       = apply_tfr_complex_baseline(S, st_baseline_time_min_ms, st_baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type);
       % abs_base_S   = abs(base_S);
       % angle_base_S = angle_S;
       
       %{
       if (j == 1)
           avg_abs_S      = abs_S;
           % avg_abs_base_S = abs_base_S;
           
           st_angle_cos = cos(angle_S);
           st_angle_sin = sin(angle_S); 
           
           % st_angle_base_cos = cos(angle_base_S);
           % st_angle_base_sin = sin(angle_base_S);
       else    
           avg_abs_S      = avg_abs_S + abs_S; 
           % avg_abs_base_S = avg_abs_base_S + abs_base_S;
           
           st_angle_cos = st_angle_cos + cos(angle_S); 
           st_angle_sin = st_angle_sin + sin(angle_S);
           
           % st_angle_base_cos = st_angle_base_cos + cos(angle_base_S); 
           % st_angle_base_sin = st_angle_base_sin + sin(angle_base_S);
       end    
       %}

    end    

    %avg_abs_S      = avg_abs_S ./ num_trials;
    avg_abs_S = squeeze(mean(abs(S_trial_vec)));
    % avg_abs_base_S = avg_abs_base_S ./ num_trials;
    
    S_trial_norm = S_trial_vec ./ abs(S_trial_vec);
    avg_angle_S = squeeze(angle(mean(S_trial_norm)));
    
    %avg_angle_S      = calc_average_phase_angle(st_angle_cos, st_angle_sin, num_trials);
    % avg_angle_base_S = calc_average_phase_angle(st_angle_base_cos, st_angle_base_sin, num_trials);
    
    [inv_S, t_inv]      = stockwell_inverse_transform(avg_abs_S, avg_angle_S, rate);
    % [inv_base_S, t_inv] = stockwell_inverse_transform(avg_abs_base_S, avg_angle_base_S, rate);

    h1_struct.data_struct.hdf1_ind_st_data(i,:)      = inv_S;
    % h1_struct.data_struct.hdf1_ind_st_base_data(i,:) = inv_base_S;

end    

% apply baseline on data

baseline_time_min_ms = h1_struct.experiment_struct.baseline_time_min_ms;
baseline_time_max_ms = h1_struct.experiment_struct.baseline_time_max_ms;
baseline_type = 1;

for j=1:n_chans
    
    x = squeeze(h1_struct.data_struct.hdf1_avg_st_data(j,:));
    y = apply_trace_baseline(x, baseline_time_min_ms, baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type);
    h1_struct.data_struct.hdf1_avg_st_data(j,1:length(y)) = y;
    
    x = squeeze(h1_struct.data_struct.hdf1_ind_st_data(j,:));
    y = apply_trace_baseline(x, baseline_time_min_ms, baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type);
    h1_struct.data_struct.hdf1_ind_st_data(j,1:length(y)) = y;
    
    % x = squeeze(h1_struct.data_struct.hdf1_ind_st_base_data(j,:));
    % y = apply_trace_baseline(x, baseline_time_min_ms, baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type);
    % h1_struct.data_struct.hdf1_ind_st_base_data(j,1:length(y)) = y;
    
end

h1_struct.st_struct.exp_name                = txt_exp;
h1_struct.st_struct.exp_case_type           = exp_case_type;
h1_struct.st_struct.the_case                = the_case;
h1_struct.st_struct.n_chans                 = n_chans;
h1_struct.st_struct.channel_list            = channel_list;
h1_struct.st_struct.max_trials              = max_trials;
h1_struct.st_struct.num_trials              = num_trials;
h1_struct.st_struct.n_samples               = size(h1_struct.data_struct.hdf1_ind_st_data,2);
h1_struct.st_struct.baseline_type           = baseline_type;
h1_struct.st_struct.st_baseline_time_min_ms = st_baseline_time_min_ms;
h1_struct.st_struct.st_baseline_time_max_ms = st_baseline_time_max_ms;
h1_struct.st_struct.hp_filter               = new_hi_pass;
h1_struct.st_struct.lp_filter               = new_lo_pass;

h1_struct.data_struct.hdf1_cnt_data   = 0;
h1_struct.data_struct.hdf1_epoch_data = 0;
	
% remove 'hdf1_prestim_epoch_data', 'hdf1_orig_data', 'hdf1_avg_data' from bloated v4-mat-files (from David, 2014-03-28)
%h1_struct.data_struct = rmfield(h1_struct.data_struct, {'hdf1_prestim_epoch_data', 'hdf1_orig_data', 'hdf1_avg_data'});

save(fname, 'h1_struct');

