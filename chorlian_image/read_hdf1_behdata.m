function [h1_struct] = read_hdf1_behdata(hdf1_filename)
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
file_struct.got_case_info = file_info.Data{6};
file_struct.got_event_info = file_info.Data{8};

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

h1_struct.case_struct = case_struct;
h1_struct.event_struct = event_struct;


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
