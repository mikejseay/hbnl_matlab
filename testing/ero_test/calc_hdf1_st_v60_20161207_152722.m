%
% function /export/home/nmanz/programs/matlab/calc_hdf1_st_v6.0.skel
%

path(path,'/export/home/kevjones/progs/hdf1progs/hdf1tools/hdf1matlab/');
path(path,'/export/home/nmanz/programs/matlab/');
path(path,'/export/home/mike/programs/matlab/');

cnthdf1_filename	= '/processed_data/cnt-h1-files/ant/256/iowa/ns32-64/ant_5_b1_30185017_256_cnt.h1';
working_directory 	= '/vol01/processed_data/mat-files-ps-v60/ns32-64/ant/ant-a/e3-n10-t100-v800';
trial_files 		= 0;	% --tf
trial_plot 			= 0;	% --tp
baseline_type 		= 1;	% -b
the_case 			= 1;	% -c
chan_list_type 		= 3;	% -e
new_hi_pass			= 0.05;	% --hi
new_lo_pass			= 55.0;	% --lo
pre_stim_time		= 500;	% -k
min_trials 			= 10;	% -n
max_trials 			= 100;	% -o
response_min_time	= 200;	% -p
response_max_time	= 0;	% -q
new_thresh_elec		= [];	% -s
threshold        	= 100;	% -t
threshold_min_time	= 0;	% -u
threshold_max_time	= 800;	% -v
st_baseline_time_min    = 0;	% -y
st_baseline_time_max    = 0;	% -z

calc_hdf1_st_v60_ms(cnthdf1_filename, working_directory, the_case, chan_list_type, baseline_type, st_baseline_time_min,st_baseline_time_max,min_trials,max_trials,pre_stim_time,threshold, threshold_min_time, threshold_max_time, new_hi_pass,new_lo_pass, new_thresh_elec, response_min_time, response_max_time, trial_plot, trial_files)
