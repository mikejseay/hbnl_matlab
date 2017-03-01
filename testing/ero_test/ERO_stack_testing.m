%% test the function /export/home/nmanz/programs/matlab/calc_hdf1_st_v60.m

% -e 3 -n 10 -t 100  -v 800

% [-e elec_list_type <1=19 default, 2=31, 3=61>]
% --> chan_list_type

% [-n out_min_trials <10 default>]  
% --> min_trials

% [-t threshold_value <-1=threshold off, 0=old value default, e.g. 100>]
% --> threshold

% [-v threshold_max_time_ms <0=old value default; real values, e.g. 600>]
% --> threshold_max_time

%cnthdf1_filename = '/processed_data/cnt-h1-files/vp3/256/suny/ns32-64/vp3_6_d1_40039067_256_cnt.h1';
%cnthdf1_filename = '/processed_data/cnt-h1-files/vp3/256/a-subjects/ns32-64/vp3_6_b2_a0000743_256_cnt.h1';
cnthdf1_filename = '/processed_data/cnt-h1-files/vp3/256/a-subjects/ns32-64/vp3_5_b1_a0001058_256_cnt.h1';

working_directory = '/active_projects/mike/test_calc';
trial_files = 0;                % --tf
trial_plot = 0;                 % --tp
baseline_type = 1;              % -b
the_case = 1;                   % -c
chan_list_type = 3;             % -e
new_hi_pass = 0.05;             % --hi
new_lo_pass = 55;               % --lo
pre_stim_time = 0;              % -k
min_trials = 10;                % -n
max_trials = 100;               % -o
response_min_time = 200;        % -p
response_max_time = 0;          % -q
%new_thresh_elec = [7,8,9,16,17,18,23,24,25]; % -s
new_thresh_elec = [];           % -s
threshold = 100;                % -t
threshold_min_time	= 0;        % -u
threshold_max_time	= 800;      % -v
st_baseline_time_min    = 0;	% -y
st_baseline_time_max    = 0;	% -z

calc_hdf1_st_v60_ms(cnthdf1_filename, working_directory, the_case, chan_list_type, ...
    baseline_type, st_baseline_time_min,st_baseline_time_max,min_trials,max_trials, ...
    pre_stim_time,threshold, threshold_min_time, threshold_max_time, new_hi_pass, ...
    new_lo_pass, new_thresh_elec, response_min_time, response_max_time, trial_plot, trial_files)

%% test the function /export/home/kevjones/progs/hdf1progs/hdf1tools/hdf1matlab/calc_chanlist_stockwell_values.m


% /export/opt/bin/solaris/extract_st_bands_$version.rb -d $dd -o new -p 1 -q 0
% -v /export/home/mort/projects/freq-for_$tw.txt -x $t1 -y $t2 -e 1 -f $tmp.3

% -d $dd = power measure <1=total default, 2=evoked, 3=induced>
    % here, changes between 1 or 2
% -o new = new output file name
    % 
% -p 1 = power units <1=power (lin) default, 2=amplitude (lin), 3=power (ln), 4=amplitude (ln)>
    % here, always 1
% -q 0 = add_baseline_values <0=no, 1=yes default>
    % here, 0
% -v /export/home/mort/projects/freq-for_$tw.txt = curve_band_type <1=all bands default, 2=single frequencies, 3=lower bands, 
				     % 4=higher bands, own freq_file>
    % here, comes from a custom text file that specifies the freqs for a
    % given timewindow. we'll use 200-500 and 4-8
% -x $t1 = min_win_time_ms <300 default>
% -y $t2 = max_win_time_ms <500 default>
% -e 1 = electrode_list_type <1=old elec_list default, 2=rows across head, 3=rows within 6 regions, own elec_list>
% -f $tmp.3 = st_mat_filelist_filename (this is a list of matfiles)

% load('test_h1_struct.mat')
% filename = '/processed_data/mat-files-v60/ns32-64/vp3/vp3-tt/e3-n10-t100-v800/a-subjects/vp3_5_b1_a0001058_256_cnt.h1.tt.st.mat';
filename = '/processed_data/mat-files-v60/ns32-64/vp3/vp3-tt/e1-n10-s9-t100-v800/suny/vp3_6_d2_40258004_256_cnt.h1.tt.st.mat';

do_baseline = 0; 										% -b (0=no default, 1=yes)
st_type 	= 1;											% -d (1=total, 2=evoked, 3=induce
% FOR V6-ALL
% channel_sort = 	1;								% -e (1=old elec_list default, 2=rows achin 6 regions, own elec_list)
% electrodes_string = '';							% -e
% electrodes_array =	[];  %		= sscanf(electrodes_string, '%f');	% -e
% FOR V6-CENTER9
channel_sort = 4;
electrodes_string = '13 16 2 10 18 5 4 19 11'; % from '/export/home/mort/projects/electrodes-for_21.txt'
electrodes_array = sscanf(electrodes_string, '%f');
% inp_files = 									% -f
file_names = {filename};
calc_type 	 = 1;										% -m (1=mean default, 2=max, 3=centroid,  6=sum)
output_text	 = 'foo.txt';        							% -o
out_type 	 = 1;										% -p (1=power (lin), 2=amplitude (lin), 3=pe (ln))
add_baseline = 0;								% -q (0=no, 1=yes default)
%frequencies_min_string 					% -v (1= all bands default, 2=single frequencieigher bands, own freq_file)	
%frequencies_max_string 				% -v
freqs_min_array 			= 4; %sscanf(frequencies_min_string, '%f');
freqs_max_array 			= 8; %sscanf(frequencies_max_string, '%f');
n_freq_files 				= min(min([length(freqs_min_array) length(freqs_max_array)]));
t_min = 300; 											% -x (300 default)
t_max = 700;											% -y (500 default)

extract_st_bands_v60_ms(do_baseline, st_type, channel_sort, electrodes_array, file_names, ...
    calc_type, output_text, out_type, add_baseline, freqs_min_array, freqs_max_array, ...
    n_freq_files, t_min, t_max)

%% putting the vals into an array

n_chans = length(S_values);
a = zeros(n_chans, 1);
for c=1:n_chans
    a(c) = S_values{1,c}{1};
end