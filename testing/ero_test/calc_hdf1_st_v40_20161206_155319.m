%
%
% function calc_hdf1_st_v4.0.skel
%
% modified by Niklas Manz (2007-08-08)
%
% adding new_hi_pass and new_lo_pass (2008-07-23)
%

path(path,'/export/home/kevjones/progs/hdf1progs/hdf1tools/hdf1matlab/');
path(path,'/export/home/nmanz/programs/matlab/');
path(path,'/export/home/mike/programs/matlab/');

cnthdf1_filename		= '/processed_data/cnt-h1-files/ans/256/c-subjects/ns16-64/ans_2_a1_c0000805_256_cnt.h1';
working_directory 	= '/vol01/processed_data/mat-files-ps-v40/ns16-64/ans/ans-f3/e1-n10-s9-t100-v800';
baseline_type 			= 1;		% -b
the_case 				= 3;		% -c
chan_list_type 		= 1;		% -e
new_hi_pass				= 0.05;		% --hi
new_lo_pass				= 55.0;		% --lo
pre_stim_time 	      = 500;		% -k
min_trials 				= 10;		% -n
max_trials 				= 100;		% -o
new_thresh_elec      = [7,8,9,16,17,18,23,24,25];       	% -s
threshold        		= 100;		% -t
threshold_min_time	= 0;		% -u
threshold_max_time	= 800;		% -v
st_baseline_time_min = 0;		% -y
st_baseline_time_max = 0;		% -z

calc_hdf1_st_v40_ms(cnthdf1_filename, working_directory, the_case, chan_list_type, baseline_type, st_baseline_time_min, st_baseline_time_max,min_trials, max_trials,pre_stim_time,threshold, threshold_min_time, threshold_max_time, new_hi_pass, new_lo_pass, new_thresh_elec)
