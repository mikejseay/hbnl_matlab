function h1_struct = calc_hdf1_st_v60(cnthdf1_filename, working_directory, the_case, chan_list_type, baseline_type, st_baseline_time_min_ms, st_baseline_time_max_ms, min_trials, max_trials, new_pre_stim_time_ms, new_threshold_value, new_threshold_min_time_ms, new_threshold_max_time_ms, new_hi_pass, new_lo_pass, new_threshold_electrodes, response_min_time, response_max_time, trial_plot, trial_files)
%
% adding of new_hi_pass and new_lo_pass variables                                           (Niklas, 2008-07-21)
% adding of new_threshold_electrodes variable                                               (Niklas, 2008-07-21)
% removed the reading of /export/age_gender_file/exp_session_subject_data.csv into h1-file  (Niklas, 2010-11-23)
% added the response time window values                                                     (Niklas, 2012-08-23)
% set the default baseline values to -180 and -50                                           (Niklas, 2012-09-13)
% added the trial_amplitude_plot option                                                     (Niklas, 2013-02-13)
% added the accepted_trial_value_files option                                               (Niklas, 2013-03-06)
% added the possibility to threshold all 61 head electrodes                                 (Niklas, 2013-05-22)
% removed the thresholding of the eye electrodes 'X' and 'Y'                                (Niklas, 2013-05-22)

today = datestr(now);

path(path,'/export/home/kevjones/progs/hdf1progs/hdf1tools/hdf1matlab/');
path(path,'/export/home/nmanz/programs/matlab/');

if (isunix == true)
    working_directory = [working_directory, '/'];
else
    working_directory = [working_directory, '\'];
end

slash = strfind(cnthdf1_filename, '/');
if (length(slash) == 0)
    slash = 0;
end

sub_id     = cnthdf1_filename(slash(length(slash))+1:length(cnthdf1_filename));
underscore = findstr(sub_id, '_');
subject_id = sub_id(underscore(3)+1:underscore(4)-1);
subject_session_run = sub_id(underscore(2)+1:underscore(2)+2);

% This is a changed program from David to run Kevin's old matlab programs with matlab2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6.1
%h1_struct = read_hdf1_data(cnthdf1_filename);
h1_struct  = read_hdf1_dataV7(cnthdf1_filename);

case_idx      = get_case_idx(h1_struct.case_struct.case_num, the_case);
exp_case_type = lower(strtrim(h1_struct.case_struct.case_type{case_idx}));

txt_exp       = cnthdf1_filename(slash(length(slash))+1:slash(length(slash))+3);
if (txt_exp == 'ans')
    if (exp_case_type == 't1' )
        exp_case_type = 'r1';
    end
    if (exp_case_type == 't2' )
        exp_case_type = 'f2';
    end
    if (exp_case_type == 't3' )
        exp_case_type = 'f3';
    end
    if (exp_case_type == 't4' )
        exp_case_type = 'f4';
    end
    if (exp_case_type == 't5' )
        exp_case_type = 'f5';
    end
    if (exp_case_type == 't6' )
        exp_case_type = 'f6';
    end
    if (exp_case_type == 't7' )
        exp_case_type = 'f7';
    end
    if (exp_case_type == 't8' )
        exp_case_type = 'f8';
    end
    if (exp_case_type == 't9' )
        exp_case_type = 'f9';
    end
end
if (txt_exp == 'aod')
    if (exp_case_type == 't' )
        exp_case_type = 'tt';
    end
end
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
    % return
end

% adding the upper limit of the response_window to the h1_struct
for i = 1:h1_struct.experiment_struct.n_cases
    if ( h1_struct.case_struct.response_window_ms(i) > 0 && response_max_time > 0 )
        h1_struct.case_struct.response_window_ms(i) = response_max_time;
    end
end

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

if (st_baseline_time_max_ms == 0 && st_baseline_time_min_ms == 0)
    st_baseline_time_min_ms = -180;
    st_baseline_time_max_ms = 0;
end
h1_struct.experiment_struct.baseline_time_min_ms = st_baseline_time_min_ms;
h1_struct.experiment_struct.baseline_time_max_ms = st_baseline_time_max_ms;

if (new_pre_stim_time_ms ~= 0)
    h1_struct.experiment_struct.pre_stim_time_ms = new_pre_stim_time_ms;
end
pre_stim_time_ms = h1_struct.experiment_struct.pre_stim_time_ms;


case_position = find(h1_struct.case_struct.case_num==the_case);
if (case_position<=size(h1_struct.case_struct.case_num,2))
    case_name       = h1_struct.case_struct.case_type{case_position};
    case_descriptor = h1_struct.case_struct.descriptor{case_position};
else
    case_name      = 'Unknown';
    case_descriptor = 'Unknown';
end

if ( h1_struct.experiment_struct.exp_name == 'ans')
    if (case_name == 'T1' )
        case_name = 'R1';
    end
    if (case_name == 'T2' )
        case_name = 'F2';
    end
    if (case_name == 'T3' )
        case_name = 'F3';
    end
    if (case_name == 'T4' )
        case_name = 'F4';
    end
    if (case_name == 'T5' )
        case_name = 'F5';
    end
    if (case_name == 'T6' )
        case_name = 'F6';
    end
    if (case_name == 'T7' )
        case_name = 'F7';
    end
    if (case_name == 'T8' )
        case_name = 'F8';
    end
    if (case_name == 'T9' )
        case_name = 'F9';
    end
end
if ( h1_struct.experiment_struct.exp_name == 'aod')
    if (case_name == 'T' )
        case_name = 'TT';
    end
end
if ( h1_struct.experiment_struct.exp_name == 'cpt')
    if (case_name == 'T' )
        case_name = 'CG';
    end
    if (case_name == 'NG' )
        case_name = 'UN';
    end
    if (case_name == 'N' )
        case_name = 'DN';
    end
end
if ( h1_struct.experiment_struct.exp_name == 'vp3')
    if (case_name == 'T' )
        case_name = 'TT';
    end
    if (case_name == 'N' )
        case_name = 'NV';
    end
end

%%%%% opening txt-file to export data of the mat-file creation process
txt_file_name = [sub_id(1:18),'case-',case_name,'.mat.v6.0.log'];
txt_exp       = cnthdf1_filename(slash(length(slash))+1:length(cnthdf1_filename)-25);

fid = fopen(txt_file_name, 'wt');
fprintf(fid, 'Information about the mat-file creation process: \n' );
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, 'Date and time:     %s \n', today);
fprintf(fid, '\n');
fprintf(fid, 'Working directory: %s \n', working_directory);
fprintf(fid, '\n');
fprintf(fid, 'Input file:        %s \n', sub_id);
fprintf(fid, '\n');
fprintf(fid, 'Subject ID:           %s \n', subject_id);
fprintf(fid, 'Experiment name:      %s \n', h1_struct.experiment_struct.exp_name);
fprintf(fid, 'Experiment session:   %s \n', subject_session_run);
fprintf(fid, 'Experiment case:      %s (%s, %s) \n', num2str(the_case), case_name, case_descriptor);
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
fprintf(fid, 'Baseline min time: %4.0f ms \n', h1_struct.experiment_struct.baseline_time_min_ms);
fprintf(fid, 'Baseline max time: %4.0f ms \n', h1_struct.experiment_struct.baseline_time_max_ms);
fprintf(fid, '\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6.2
h1_struct = create_cnthdf1_baseline_epoch_data_v6(h1_struct, txt_file_name, new_hi_pass, new_lo_pass, new_threshold_electrodes, the_case, trial_plot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6.3
h1_struct = determine_correct_and_accepted_trials_v6(h1_struct, txt_file_name, case_name, response_min_time, response_max_time);

fclose(fid);

trial_file_name = [sub_id(1:17),'_event+trial-info.v6.0.log'];
trialfile = fopen(trial_file_name, 'wt');
fprintf(trialfile, 'Information about the event and trial structure: \n');
fprintf(trialfile, '\n');
fprintf(trialfile, '\n');
fprintf(trialfile, 'Date:               %s\n', date);
fprintf(trialfile, '\n');
fprintf(trialfile, 'Working directroy: %s\n', working_directory);
fprintf(trialfile, 'Input file:        %s\n', sub_id);
fprintf(trialfile, '\n');
fprintf(trialfile, 'Subject_ID: %s\n', subject_id);
fprintf(trialfile, 'Experiment_name: %s\n', h1_struct.experiment_struct.exp_name);
fprintf(trialfile, 'Experiment_session: %s\n', subject_session_run);
fprintf(trialfile, '\n');
fprintf(trialfile, 'Threshold_min_time: %5.1f\n', h1_struct.experiment_struct.threshold_time_min_ms);
fprintf(trialfile, 'Threshold_max_time: %5.1f\n', h1_struct.experiment_struct.threshold_time_max_ms);
fprintf(trialfile, 'Threshold_value: %3.0f\n',    h1_struct.experiment_struct.threshold_value);
fprintf(trialfile, 'lp_filter: %5.2f\n', new_lo_pass);
fprintf(trialfile, 'hp_filter: %3.2f\n', new_hi_pass);
fprintf(trialfile, 'Baseline_min_time: %3.0f\n', h1_struct.experiment_struct.baseline_time_min_ms);
fprintf(trialfile, 'Baseline_max_time: %3.0f\n', h1_struct.experiment_struct.baseline_time_max_ms);
fprintf(trialfile, 'Response_window_min_time: %4.0f\n', response_min_time);
fprintf(trialfile, 'Response_window_max_time: %4.0f\n', h1_struct.case_struct.response_window_ms(case_position));
fprintf(trialfile, '\n');
fprintf(trialfile, 'events: %3.0f\n',                                                 h1_struct.experiment_struct.n_events);
fprintf(trialfile, 'event_numbers: %s\n', num2str(h1_struct.event_struct.event_num( 1:h1_struct.experiment_struct.n_events)));
fprintf(trialfile, 'event_id: %s\n',      num2str(h1_struct.event_struct.event_id(  1:h1_struct.experiment_struct.n_events)));
fprintf(trialfile, 'keypad_id: %s\n',     num2str( h1_struct.event_struct.keypad_id(1:h1_struct.experiment_struct.n_events)));
fprintf(trialfile, '\n');
fprintf(trialfile, 'trials: %3.0f\n',                                                              h1_struct.experiment_struct.n_trials);
fprintf(trialfile, 'trial_number: %s\n',     num2str(h1_struct.trial_struct.trial_num(           1:h1_struct.experiment_struct.n_trials)));
fprintf(trialfile, 'case_number: %s\n',      num2str(h1_struct.trial_struct.case_num(            1:h1_struct.experiment_struct.n_trials)));
fprintf(trialfile, 'response_id: %s\n',      num2str(h1_struct.trial_struct.response_id(         1:h1_struct.experiment_struct.n_trials)));
fprintf(trialfile, 'stim_id: %s\n',          num2str(h1_struct.trial_struct.stim_id(             1:h1_struct.experiment_struct.n_trials)));
fprintf(trialfile, 'correct: %s\n',          num2str(h1_struct.trial_struct.correct(             1:h1_struct.experiment_struct.n_trials)));
fprintf(trialfile, 'omitted: %s\n',          num2str(h1_struct.trial_struct.omitted(             1:h1_struct.experiment_struct.n_trials)));
fprintf(trialfile, 'artifact: %s\n',         num2str(h1_struct.trial_struct.artf_present(        1:h1_struct.experiment_struct.n_trials)));
fprintf(trialfile, 'accepted: %s\n',         num2str(h1_struct.trial_struct.accepted(            1:h1_struct.experiment_struct.n_trials)));
fprintf(trialfile, 'max_window_value:');
fprintf(trialfile, '%10.0f', h1_struct.trial_struct.max_thresh_win_value(1:h1_struct.experiment_struct.n_trials));
fprintf(trialfile, '\n');
fprintf(trialfile, 'response_latency:');
fprintf(trialfile, '%5.0f', h1_struct.trial_struct.response_latency(1:h1_struct.experiment_struct.n_trials));
fprintf(trialfile, '\n');
fclose(trialfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6.4
h1_struct = create_epoch_average_data_v6(h1_struct);

clear h1_struct.data_struct.hdf1_ind_st_data;

n_trials  = size(h1_struct.data_struct.hdf1_epoch_data,1);
n_samples = size(h1_struct.data_struct.hdf1_epoch_data,3);
n_cases   = h1_struct.experiment_struct.n_cases;
rate      = h1_struct.experiment_struct.rate;

fid = fopen(txt_file_name, 'a');
fprintf(fid, '\n');
fprintf(fid, 'Case numbers               :');
for i = 1:n_cases
    fprintf(fid, '%3.0f ', h1_struct.case_struct.case_num(i));
end
fprintf(fid, '\n');
fprintf(fid, 'Case names                 :');
for i = 1:n_cases
    if ( h1_struct.experiment_struct.exp_name == 'ans')
        if (h1_struct.case_struct.case_type{i} == 'T1' )
            h1_struct.case_struct.case_type{i} = 'R1';
        end
        if (h1_struct.case_struct.case_type{i} == 'T2' )
            h1_struct.case_struct.case_type{i} = 'F2';
        end
        if (h1_struct.case_struct.case_type{i} == 'T4' )
            h1_struct.case_struct.case_type{i} = 'F4';
        end
        if (h1_struct.case_struct.case_type{i} == 'T5' )
            h1_struct.case_struct.case_type{i} = 'F5';
        end
        if (h1_struct.case_struct.case_type{i} == 'T6' )
            h1_struct.case_struct.case_type{i} = 'F6';
        end
        if (h1_struct.case_struct.case_type{i} == 'T7' )
            h1_struct.case_struct.case_type{i} = 'F7';
        end
        if (h1_struct.case_struct.case_type{i} == 'T8' )
            h1_struct.case_struct.case_type{i} = 'F8';
        end		
        if (h1_struct.case_struct.case_type{i} == 'T9' )
            h1_struct.case_struct.case_type{i} = 'F9';
        end		
    end
    if ( h1_struct.experiment_struct.exp_name == 'aod')
        if (h1_struct.case_struct.case_type{i} == 'T' )
            h1_struct.case_struct.case_type{i} = 'TT';
        end
    end
    if ( h1_struct.experiment_struct.exp_name == 'cpt')
        if (h1_struct.case_struct.case_type{i} == 'T' )
            h1_struct.case_struct.case_type{i} = 'CG';
        end
        if (h1_struct.case_struct.case_type{i} == 'NG' )
            h1_struct.case_struct.case_type{i} = 'UN';
        end
        if (h1_struct.case_struct.case_type{i} == 'N' )
            h1_struct.case_struct.case_type{i} = 'DN';
        end
    end
    if ( h1_struct.experiment_struct.exp_name == 'vp3')
        if (h1_struct.case_struct.case_type{i} == 'T' )
            h1_struct.case_struct.case_type{i} = 'TT';
        end
        if (h1_struct.case_struct.case_type{i} == 'N' )
            h1_struct.case_struct.case_type{i} = 'NV';
        end
    end

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
fprintf(fid, 'Number of trials (artifact):');
for i = 1:n_cases
    fprintf(fid, '%3.0f ', h1_struct.case_struct.n_artfs(i));
end
fprintf(fid, '\n');
fprintf(fid, 'Number of trials (accepted):');
for i = 1:n_cases
    fprintf(fid, '%3.0f ', h1_struct.case_struct.n_trials_accepted(i));
end
fprintf(fid, '\n');
fprintf(fid, '\n');

% checking if number of accepted trials is larger than number of responses (if response required) and print note into log file
for i = 1:n_cases
    if ( h1_struct.case_struct.response_id(i) > 0 )
        if ( h1_struct.case_struct.n_responses(i)-h1_struct.case_struct.n_resp_errors(i) > 0)
            if (h1_struct.case_struct.n_trials_accepted(i) > h1_struct.case_struct.n_responses(i))
                fprintf(fid, 'Number of accepted trials is larger than number of responses (if required) for the %s case!!\n',  h1_struct.case_struct.case_type{i});
                fprintf(fid, '\n');
                fprintf(fid, '\n');
            end
        end
    end
end

% checking if number of accepted trials is larger than number of responses (if response required), then return
if ( h1_struct.case_struct.n_responses(the_case)-h1_struct.case_struct.n_resp_errors(the_case) > 0)
    if (h1_struct.case_struct.n_trials_accepted(the_case) > h1_struct.case_struct.n_responses(the_case))
        return
    end
end

if (case_idx == 0)
    return
end

case_correct_idx  = get_correct_case_trial_idxs(h1_struct, the_case);  %%%%% 6.5
case_omitted_idx  = get_omitted_case_trial_idxs(h1_struct, the_case);  %%%%% 6.6
case_artifact_idx = get_artifact_case_trial_idxs(h1_struct, the_case); %%%%% 6.7
case_accepted_idx = get_accepted_case_trial_idxs(h1_struct, the_case); %%%%% 6.8

% adding the part to print the accepted trial values into a separate file
if (trial_files == 1 )
    txt_exp         = cnthdf1_filename(slash(length(slash))+1:length(cnthdf1_filename)-25);
    for i = 1:length(case_accepted_idx)
        trial_num       = case_accepted_idx(i);
        case_num        = h1_struct.trial_struct.case_num(trial_num);
        trial_file_name = [sub_id(1:18),'trial-',num2str(trial_num),'_case-',num2str(case_num),'.txt'];
        fid2 = fopen(trial_file_name, 'wt');
        fprintf(fid2, 'sample# ');
        for k=1:length(h1_struct.run_struct.channel_label)
            fprintf(fid2, '%s ', h1_struct.run_struct.channel_label{k});
        end
        fprintf(fid2, '\n');
        for j=1:length(h1_struct.data_struct.hdf1_epoch_data(1,1,:))
            fprintf(fid2, '%s %s \n', num2str(j), num2str(h1_struct.data_struct.hdf1_epoch_data(i,:,j)) );
        end
        fclose(fid2);
    end
end

fprintf(fid, 'Trial information for case %s (%s):\n', num2str(the_case), case_name);
fprintf(fid, 'Number of correct trials       : %s \n',  num2str(length(case_correct_idx)));
fprintf(fid, 'Correct trials                 : %s \n',  num2str(case_correct_idx));
fprintf(fid, 'Number of omitted trials       : %s \n',  num2str(length(case_omitted_idx)));
if ( length(case_omitted_idx) == 0 )
    fprintf(fid, 'Omitted trials                 : none \n');
elseif ( length(case_omitted_idx) == 1 )
    fprintf(fid, 'Omitted trial                  : %s \n',  num2str(case_omitted_idx));
else
    fprintf(fid, 'Omitted trials                 : %s \n',  num2str(case_omitted_idx));
end
fprintf(fid, 'Number of artifact trials      : %s \n',  num2str(length(case_artifact_idx)));
fprintf(fid, 'Artifact trials                : %s \n',  num2str(case_artifact_idx));
fprintf(fid, 'Number of accepted trials      : %s \n',  num2str(length(case_accepted_idx)));
fprintf(fid, 'Accepted trials                : %s \n',  num2str(case_accepted_idx));
case_accepted_idx = case_accepted_idx(randperm(length(case_accepted_idx)));

num_trials = min([length(case_accepted_idx) max_trials]);
fprintf(fid, '\n');
fprintf(fid, 'Defined number of max trials   : %s \n',  num2str(max_trials));

% Check if number of accepted trials is not equal to 0 or smaller than the minimal number of trials.
if (length(case_accepted_idx) == 0)
    fprintf(fid, '\n');
    fprintf(fid, '!!!!! No accepted trials for ERO calculation !!!!! \n');
    return
end

if (length(case_accepted_idx) < min_trials)
    fprintf(fid, '\n');
    fprintf(fid, '!!!!! Not enough accepted trials for ERO calculation !!!!! \n');
    return
end

fprintf(fid, 'Number of used ERO trials      : %s \n',  num2str(length(case_accepted_idx(1:num_trials))));
fprintf(fid, 'Used trials for ERO calculation: %s \n',  num2str(case_accepted_idx(1:num_trials)));
fclose(fid);

[channel_list, n_chans] = get_channel_list(chan_list_type); %%%%% 6.9
channel_indices = get_channel_indices(h1_struct.run_struct.channel_label, channel_list); %%%%% 6.10

Wn_hi = new_hi_pass / (rate/2);
%[bhi,ahi] = cheby1(3,0.5,Wn_hi,'high');
[bhi,ahi] = cheby1(5,0.043,Wn_hi,'high');


for i=1:n_chans
    channel_idx = channel_indices(i);
    channel_name = get_channel_name(h1_struct.run_struct.channel_label, channel_idx); %%%%% 6.11

    sprintf('%s', [num2str(i), ': ', channel_name{:}]);

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
    mean_data = filtfilt(bhi,ahi,mean_data);
    [Savg,t,f] = stockwell_transform(mean_data, 0, length(mean_data)/2, 1/rate); %%%%% 6.12
    abs_Savg   = abs(Savg);
    angle_Savg = angle(Savg);

    for j=1:num_trials
        trial_vector = trial_data(j,:);
        trial_vector = filtfilt(bhi,ahi,trial_vector);
        [S,t,f] = stockwell_transform(trial_vector, 0, length(trial_vector)/2, 1/rate); %%%%% 6.12

        abs_S   = abs(S);
        angle_S = angle(S);

        base_S       = apply_tfr_complex_baseline(S, st_baseline_time_min_ms, st_baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type); %%%%% 6.13
        abs_base_S   = abs(base_S);
        angle_base_S = angle_S;

        if (j == 1)
            avg_abs_S      = abs_S;
            avg_abs_base_S = abs_base_S;

            st_angle_cos = cos(angle_S);
            st_angle_sin = sin(angle_S);

            st_angle_base_cos = cos(angle_base_S);
            st_angle_base_sin = sin(angle_base_S);
        else
            avg_abs_S      = avg_abs_S + abs_S;
            avg_abs_base_S = avg_abs_base_S + abs_base_S;

            st_angle_cos = st_angle_cos + cos(angle_S);
            st_angle_sin = st_angle_sin + sin(angle_S);

            st_angle_base_cos = st_angle_base_cos + cos(angle_base_S);
            st_angle_base_sin = st_angle_base_sin + sin(angle_base_S);
        end
    end

    avg_abs_S      = avg_abs_S ./ num_trials;
    avg_abs_base_S = avg_abs_base_S ./ num_trials;

    avg_angle_S      = calc_average_phase_angle(st_angle_cos, st_angle_sin, num_trials); %%%%% 6.14
    avg_angle_base_S = calc_average_phase_angle(st_angle_base_cos, st_angle_base_sin, num_trials); %%%%% 6.14

    [inv_S, t_inv]      = stockwell_inverse_transform(avg_abs_S, avg_angle_S, rate);
    [inv_base_S, t_inv] = stockwell_inverse_transform(avg_abs_base_S, avg_angle_base_S, rate);

    h1_struct.data_struct.hdf1_ind_st_data(i,:)      = inv_S;
    h1_struct.data_struct.hdf1_ind_st_base_data(i,:) = inv_base_S;
end

% apply baseline on data

baseline_time_min_ms = h1_struct.experiment_struct.baseline_time_min_ms;
baseline_time_max_ms = h1_struct.experiment_struct.baseline_time_max_ms;
baseline_type = 1;

if (abs(baseline_time_min_ms) > pre_stim_time_ms )
    baseline_time_min_ms = -pre_stim_time_ms;
end


for j=1:n_chans

    x = squeeze(h1_struct.data_struct.hdf1_avg_st_data(j,:));
    y = apply_trace_baseline(x, baseline_time_min_ms, baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type);
    h1_struct.data_struct.hdf1_avg_st_data(j,1:length(y)) = y;

    x = squeeze(h1_struct.data_struct.hdf1_ind_st_data(j,:));
    y = apply_trace_baseline(x, baseline_time_min_ms, baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type);
    h1_struct.data_struct.hdf1_ind_st_data(j,1:length(y)) = y;

    x = squeeze(h1_struct.data_struct.hdf1_ind_st_base_data(j,:));
    y = apply_trace_baseline(x, baseline_time_min_ms, baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type);
    h1_struct.data_struct.hdf1_ind_st_base_data(j,1:length(y)) = y;

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

%save(fname, 'h1_struct');

