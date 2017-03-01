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

cnthdf1_filename		= '/vol01/processed_data/cnt-h1-files/vp3/256/suny/ns32-64/vp3_5_b1_40056125_256_cnt.h1';
working_directory 	= '/export/home/mike/ero_test';
baseline_type 			= 1;		% -b
the_case 				= 1;		% -c
chan_list_type 		= 1;		% -e
new_hi_pass				= 0.05;		% --hi
new_lo_pass				= 55.0;		% --lo
pre_stim_time 	      = 500;		% -k
min_trials 				= 10;		% -n
max_trials 				= 100;		% -o
new_thresh_elec      = [2,4,5,10,11,13,16,18,19];       	% -s
threshold        		= 100;		% -t
threshold_min_time	= 0;		% -u
threshold_max_time	= 800;		% -v
st_baseline_time_min = 0;		% -y
st_baseline_time_max = 0;		% -z

calc_hdf1_st_v40_ms(cnthdf1_filename, working_directory, the_case, chan_list_type, baseline_type, st_baseline_time_min, st_baseline_time_max,min_trials, max_trials,pre_stim_time,threshold, threshold_min_time, threshold_max_time, new_hi_pass, new_lo_pass, new_thresh_elec)
