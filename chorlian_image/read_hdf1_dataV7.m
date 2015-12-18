function [h1_struct] = read_hdf1_dataV7(hdf1_filename)
% see important note at end of file
% function [h1_struct] = read_hdf1_dataV7(hdf1_filename)
% for use with 
hdf1_int_data_path = [];
v = version('-release');
if length(v) < 4
	h1_struct = read_hdf1_data(hdf1_filename);
	return
end

hdf1_info = hdf5info(hdf1_filename, 'V71Dimensions', true);

% read relevant file info

file_info = hdf5read(hdf1_info.GroupHierarchy.Groups(1).Datasets(1), 'V71Dimensions', true);
file_struct.data_type = file_info.Data{1};
file_struct.original_file = file_info.Data{2};
file_struct.got_subject_info = file_info.Data{3};
file_struct.got_exp_info = file_info.Data{4};
file_struct.got_run_info = file_info.Data{5};
file_struct.got_case_info = file_info.Data{6};
file_struct.got_trial_info = file_info.Data{7};
file_struct.got_event_info = file_info.Data{8};
file_struct.got_transforms_info = file_info.Data{9};
file_struct.got_data = file_info.Data{10};
file_struct.n_dimensions = file_info.Data{13};
file_struct.n_subjects = file_info.Data{14};
file_struct.file_name = file_info.Data{16}.Data.Data;

% read in subject info

if (file_struct.got_subject_info == 1)
 subject_info = hdf5read(hdf1_info.GroupHierarchy.Groups(1).Groups(3).Datasets(1), 'V71Dimensions', true);
 num_subjects = size(subject_info,2);
 for i = 1:num_subjects
 exp_version(i) = subject_info(i).Data(1);
 run_number(i) = subject_info(i).Data(2);
 age = subject_info(i).Data(3);
 subject_id{i} = subject_info(i).Data{4}.Data.Data;
 dob{i} = subject_info(i).Data{5}.Data.Data;
 gender{i} = subject_info(i).Data{6}.Data.Data;
 handedness{i} = subject_info(i).Data{7}.Data.Data;
 session_code{i} = subject_info(i).Data{8}.Data.Data;
 end
 
 subject_struct.exp_version = exp_version;
 subject_struct.run_number = run_number;
 subject_struct.age = age;
 subject_struct.subject_id = subject_id;
 subject_struct.dob = dob;
 subject_struct.gender = gender;
 subject_struct.handedness = handedness;
 subject_struct.session_code = session_code;
else
 subject_struct = 0;
end


 
% read relevant experiment info

if (file_struct.got_exp_info == 1)
 experiment_info = hdf5read(hdf1_info.GroupHierarchy.Groups(1).Groups(1).Datasets(1),  'V71Dimensions',  true);
 experiment_struct.n_cases = experiment_info.Data{2};
 experiment_struct.n_trials = experiment_info.Data{3};
 experiment_struct.n_events = experiment_info.Data{4};
 experiment_struct.n_chans = experiment_info.Data{5};
 experiment_struct.n_samples = experiment_info.Data{6};
 experiment_struct.n_dim1 = experiment_info.Data{7};
 experiment_struct.n_lines = experiment_info.Data{8};
 experiment_struct.n_transforms = experiment_info.Data{9};
 experiment_struct.exp_version = experiment_info.Data{10};
 experiment_struct.rate = double(experiment_info.Data{11});
 experiment_struct.apply_thresholding = double(experiment_info.Data{12});
 experiment_struct.avg_trigisi_ms = double(experiment_info.Data{13});
 experiment_struct.pre_stim_time_ms = experiment_info.Data{14};
 experiment_struct.post_stim_time_ms = experiment_info.Data{15};
 experiment_struct.baseline_time_min_ms = experiment_info.Data{16};
 experiment_struct.baseline_time_max_ms = experiment_info.Data{17};
 experiment_struct.samp_interval = experiment_info.Data{18};
 experiment_struct.uv_per_unit = experiment_info.Data{19};
 experiment_struct.threshold_time_min_ms = experiment_info.Data{20};
 experiment_struct.threshold_time_max_ms = experiment_info.Data{21};
 experiment_struct.threshold_value = experiment_info.Data{22};
  experiment_struct.exp_name = experiment_info.Data{23}.Data.Data;
else
 experiment_struct = [];
end 

% read in relevant run data 

if (file_struct.got_run_info == 1)
 run_info = hdf5read(hdf1_info.GroupHierarchy.Groups(1).Groups(2).Datasets(1), 'V71Dimensions', true);
 run_struct.run_number = run_info.Data{1};
 run_struct.dc_baseline = run_info.Data{2}.Data;
 run_struct.sensitivity = run_info.Data{3}.Data;
 run_struct.calibration = run_info.Data{4}.Data;
 run_struct.run_date_time = run_info.Data{6}.Data.Data;
 
 num_chans = size(run_info.Data{7}.Data,2);
 for i = 1:num_chans
 channel_label{i} = run_info.Data{7}.Data(i).Data;
 end 
 run_struct.channel_label = channel_label;
 
 run_struct.session_code = run_info.Data{8}.Data.Data;
 
 if (file_struct.original_file == 1) 
 scale_value = zeros(1,num_chans) + experiment_struct.uv_per_unit;
 run_struct.scale_value = scale_value;
 else
 scale_value = run_struct.sensitivity.*(run_struct.calibration./204.8);
 run_struct.scale_value = scale_value;
 experiment_struct.uv_per_unit = scale_value; 
 end 
 
else
 run_struct = [];
end 

% read in data 

if (file_struct.got_data == 1)
 if (file_struct.data_type == 0)
 data_struct.hdf1_data_dims = hdf1_info.GroupHierarchy.Groups(2).Groups(1).Datasets(1).Dims;
 hdf1_int_data = hdf5read(hdf1_info.GroupHierarchy.Groups(2).Groups(1).Datasets(1), 'V71Dimensions', true);
 data_struct.hdf1_int_data_class = class(hdf1_int_data);
 data_struct.hdf1_int_data_path = hdf1_info.GroupHierarchy.Groups(2).Groups(1).Datasets(1).Name;
 data_struct.hdf1_cnt_data = double(hdf1_int_data) * double(experiment_struct.uv_per_unit(1));
 data_struct.hdf1_epoch_data = 0;
 data_struct.hdf1_avg_data = 0;
 else
 data_struct.hdf1_data_dims = hdf1_info.GroupHierarchy.Groups(2).Datasets(1).Dims;
 data_struct.hdf1_cnt_data = 0;
 data_struct.hdf1_epoch_data = 0;
 hdf1_avg_data = hdf5read(hdf1_info.GroupHierarchy.Groups(2).Datasets(1), 'V71Dimensions', true);
 
 % re-organise avg data
 
 for the_case=1:size(hdf1_avg_data,3) 
 for channel=1:size(hdf1_avg_data,1) 
 for sample=1:size(hdf1_avg_data,2)
 data_struct.hdf1_avg_data(the_case,channel,sample) = double(hdf1_avg_data(channel,sample,the_case) * experiment_struct.uv_per_unit(1));
 end 
 end
 end 
 
 end
else
 data_struct = [];
end

% read in relevant case info

if (file_struct.got_case_info == 1)
 case_info = hdf5read(hdf1_info.GroupHierarchy.Groups(1).Groups(2).Groups(1).Datasets(1), 'V71Dimensions', true);
 num_cases = size(case_info,2);
 for i = 1:num_cases
 case_nums(i) = case_info(i).Data{1}.Data;
 n_trials(i) = case_info(i).Data{2}.Data;
 n_trials_accepted(i) = case_info(i).Data{3}.Data;
 stim_id(i) = case_info(i).Data{4}.Data;
 response_id(i) = case_info(i).Data{5}.Data;
 n_artfs(i) = case_info(i).Data{6}.Data;
 n_responses(i) = case_info(i).Data{7}.Data;
 n_resp_errors(i) = case_info(i).Data{8}.Data;
 n_resp_omitted(i) = case_info(i).Data{9}.Data;
 response_window_ms(i) = case_info(i).Data{10}.Data;
 mean_resp_time(i) = case_info(i).Data{11}.Data;
 mean_resp_accuracy(i) = case_info(i).Data{12}.Data;
 case_type{i} = case_info(i).Data{14}.Data.Data;
 descriptor{i} = case_info(i).Data{15}.Data.Data;
 end
 case_struct.case_nums = case_nums;
 case_struct.case_num = case_nums;
 case_struct.n_trials = n_trials;
 case_struct.n_trials_accepted = n_trials_accepted;
 case_struct.stim_id = stim_id;
 case_struct.response_id = response_id;
 case_struct.n_artfs = n_artfs;
 case_struct.n_responses = n_responses;
 case_struct.n_resp_errors = n_resp_errors;
 case_struct.n_resp_omitted = n_resp_omitted;
 case_struct.response_window_ms = response_window_ms;
 case_struct.mean_resp_time = mean_resp_time;
 case_struct.mean_resp_accuracy = mean_resp_accuracy;
 case_struct.case_type = case_type;
 case_struct.descriptor = descriptor;
else
 case_struct = [];
end 

% read in relevant event info

if (file_struct.got_event_info == 1)
 event_info = hdf5read(hdf1_info.GroupHierarchy.Groups(1).Groups(2).Groups(2).Datasets(1), 'V71Dimensions', true);
 num_events = size(event_info,2);
 for i = 1:num_events
 event_num(i) = event_info(i).Data{1}.Data;
 event_id(i) = event_info(i).Data{2}.Data;
 keypad_id(i) = event_info(i).Data{3}.Data;
 keyboard_id(i) = event_info(i).Data{4}.Data;
 mouse_id(i) = event_info(i).Data{5}.Data;
 file_offset(i) = event_info(i).Data{6}.Data;
 event_time_offset(i) = event_info(i).Data{7}.Data;
 end 
 event_struct.event_num = event_num;
 event_struct.event_id = event_id;
 event_struct.keypad_id = keypad_id;
 event_struct.keyboard_id = keyboard_id;
 event_struct.mouse_id = mouse_id;
 event_struct.file_offset = file_offset;
 event_struct.event_time_offset = event_time_offset;
else
 event_struct = [];
end 

% read in relevant trial info

if (file_struct.got_trial_info == 1)
 trial_info = hdf5read(hdf1_info.GroupHierarchy.Groups(1).Groups(2).Groups(3).Datasets(1), 'V71Dimensions', true);
 num_trials = size(trial_info,2);
 for i = 1:num_trials
 trial_num(i) = trial_info(i).Data{1}.Data;
 case_num(i) = trial_info(i).Data{2}.Data;
 response_id(i) = trial_info(i).Data{3}.Data;
 stim_id(i) = trial_info(i).Data{4}.Data;
 correct(i) = trial_info(i).Data{5}.Data;
 omitted(i) = trial_info(i).Data{6}.Data;
 artf_present(i) = trial_info(i).Data{7}.Data;
 accepted(i) = trial_info(i).Data{8}.Data;
 max_thresh_win_value(i) = trial_info(i).Data{9}.Data;
 threshold_value(i) = trial_info(i).Data{10}.Data;
 response_latency(i) = trial_info(i).Data{11}.Data;
 time_offset(i) = trial_info(i).Data{12}.Data;
 end 
 trial_struct.trial_num = trial_num;
 trial_struct.case_num = case_num;
 trial_struct.response_id = response_id;
 trial_struct.stim_id= stim_id;
 trial_struct.correct = correct;
 trial_struct.omitted= omitted;
 trial_struct.artf_present = artf_present;
 trial_struct.accepted = accepted;
 trial_struct.max_thresh_win_value =max_thresh_win_value;
 trial_struct.threshold_value= threshold_value;
 trial_struct.response_latency =response_latency;
 trial_struct.time_offset = time_offset;
else
 trial_struct = [];
end 

% read in relevant trial info

if (file_struct.got_transforms_info == 1)
 transforms_info = hdf5read(hdf1_info.GroupHierarchy.Groups(1).Groups(4).Datasets(1), 'V71Dimensions', true);
 transforms_struct.hi_pass_filter = transforms_info.Data{1}.Data;
 transforms_struct.lo_pass_filter = transforms_info.Data{2}.Data;
else
 transforms_struct = [];
end 

h1_struct.file_struct = file_struct;
h1_struct.subject_struct = subject_struct;
h1_struct.data_struct = data_struct;
h1_struct.experiment_struct = experiment_struct;
h1_struct.run_struct = run_struct;
h1_struct.case_struct = case_struct;
h1_struct.event_struct = event_struct;
h1_struct.trial_struct = trial_struct;
h1_struct.transforms_struct = transforms_struct;


% >> hinfo = hdf5info('myfile.h5');
%      >> data = hinfo.GroupHierarchy.Datasets(1);
%      >> hdf5path = hinfo.GroupHierarchy.Datasets(1).Name;
% 
% The variable 'data' will contain the data necessary for your data analysis, and 
% 'hdf5path' will contain the field information to write that data into an hdf5path. 
% 
% The syntax to do this with the HDF5WRITE function would then be as follows: 
% 
%      >> hdf5write('myresults.h5',hdf5path,data);
