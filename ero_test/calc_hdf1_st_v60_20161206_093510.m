%
% function /export/home/nmanz/programs/matlab/calc_hdf1_st_v6.0.skel
%

path(path,'/export/home/kevjones/progs/hdf1progs/hdf1tools/hdf1matlab/');
path(path,'/export/home/nmanz/programs/matlab/');
path(path,'/export/home/mike/programs/matlab/');

cnthdf1_filename	= '/processed_data/cnt-h1-files/vp3/256/a-subjects/mc16-64/vp3_2_a1_a0000479_256_cnt.h1';
working_directory 	= '/vol01/processed_data/mat-files-ps-v60/mc16-64/vp3/vp3-nt/e1-n10-s9-t100-v800';
trial_files 		= 0;	% --tf
trial_plot 			= 0;	% --tp
baseline_type 		= 1;	% -b
the_case 			= 2;	% -c
chan_list_type 		= 1;	% -e
new_hi_pass			= 0.05;	% --hi
new_lo_pass			= 55.0;	% --lo
pre_stim_time		= 600;	% -k
min_trials 			= 10;	% -n
max_trials 			= 100;	% -o
response_min_time	= 200;	% -p
response_max_time	= 0;	% -q
new_thresh_elec		= [7,8,9,16,17,18,23,24,25];	% -s
threshold        	= 100;	% -t
threshold_min_time	= 0;	% -u
threshold_max_time	= 800;	% -v
st_baseline_time_min    = 0;	% -y
st_baseline_time_max    = 0;	% -z

calc_hdf1_st_v60_ms(cnthdf1_filename, working_directory, the_case, chan_list_type, baseline_type, st_baseline_time_min,st_baseline_time_max,min_trials,max_trials,pre_stim_time,threshold, threshold_min_time, threshold_max_time, new_hi_pass,new_lo_pass, new_thresh_elec, response_min_time, response_max_time, trial_plot, trial_files);
% h1_struct_nm = calc_hdf1_st_v60(cnthdf1_filename, working_directory, the_case, chan_list_type, baseline_type, st_baseline_time_min,st_baseline_time_max,min_trials,max_trials,pre_stim_time,threshold, threshold_min_time, threshold_max_time, new_hi_pass,new_lo_pass, new_thresh_elec, response_min_time, response_max_time, trial_plot, trial_files);
