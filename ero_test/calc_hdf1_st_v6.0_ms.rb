#! /usr/local/bin/ruby
#
# /opt/bin/solaris/calc_hdf1_st_v6.0.rb
#
# Niklas Manz (???)
#

require 'getoptlong'
require 'ftools'
require '/export/home/nmanz/programs/ruby/rubytools_2a'

usage = "#{$0}\n
	This program will check if a mat-file for an ID already exists.
	It compares ONLY the file name not the internal parameter!
	If the mat-file does not exist in the folder you are in, it will call calc_hdf1_st_v4\n
	[--tf accepted_trial_value_files <0=no default, 1=yes>] 
	[--tp trial_amplitude_plot <0=no default, 1=yes>] 
	[-b norm_baseline_type_to_use <1=mean default, 2=rms, 3=geomean, 4=harmmean>] 
	[-c case_number <e.g., 1=tt, 2=nt, 3=nv>] 
	[-e elec_list_type <1=19 default, 2=31, 3=61>] 
	[-g use_grid <0=no default, 1=yes>] 
	[--hi hp_filter <0.05 default>] 
	[--lo lp_filter <55.0 default>]
	[-k pre_stim_time_ms <positive value!>] 
	[-n out_min_trials <10 default>] 
	[-o out_max_trials <100 default>] 
	[-p response_window_min_time_ms <200 default>]
	[-q response_window_max_time_ms <0=old value default>]
	[-s threshold_electrodes <first 31 default, 61, own comma separated electrode list>] 
	[-t threshold_value <-1=threshold off, 0=old value default, e.g. 100>] 
	[-u threshold_min_time_ms <0=old value default; real values, e.g. -100>] 
	[-v threshold_max_time_ms <0=old value default; real values, e.g. 600>] 
	[-y st_baseline_time_min_ms <real values, i.e., negative if pre stimulus>] 
	[-z st_baseline_time_max_ms <real values; i.e., negative if pre stimulus>]\n
	[-f filelist_filename]\n"

# specify the options we accept and initialize
# the option parser

opts = GetoptLong.new(	[ "--trial_amplitude_plot",		"--tp",GetoptLong::REQUIRED_ARGUMENT ],
			[ "--accepted_trial_value_files",	"--tf",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--norm_baseline_type_to_use",	"-b",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--case_number", 			"-c",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--elec_list_type", 			"-e",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--filelist_filename", 		"-f",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--use_grid", 			"-g",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--hp_filter",	 		"--hi",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--lp_filter",	 		"--lo",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--prepended_id_list", 		"-i",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--pre_stim_time_ms",        		"-k",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--out_min_trials", 	       		"-n",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--out_max_trials", 	       		"-o",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--response_window_min_time_ms",	"-p",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--response_window_max_time_ms",	"-q",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--threshold_electrodes",    		"-s",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--threshold_value", 	       		"-t",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--threshold_min_time_ms",   		"-u",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--threshold_max_time_ms",   		"-v",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--st_baseline_time_min_ms", 		"-y",	GetoptLong::REQUIRED_ARGUMENT ],
			[ "--st_baseline_time_max_ms", 		"-z",	GetoptLong::REQUIRED_ARGUMENT ]
)

# process the parsed options

filelist_filename = ''
prepended_id_list = 0

trial_amplitude_plot		= 0
accepted_trial_value_files	= 0
norm_baseline_type_to_use 	= 1
case_number 			= 1
elec_list_type 			= 1
use_grid 			= 0
hp_filter	 		= 0.05
lp_filter			= 55.0
pre_stim_time_ms 		= 0
out_min_trials 			= 10
out_max_trials 			= 100
response_window_min_time_ms	= 200
response_window_max_time_ms	= 0
threshold_value 		= 75
threshold_electrodes 		= ""
threshold_min_time_ms 		= 0 
threshold_max_time_ms 		= 0 
st_baseline_time_min_ms 	= 0
st_baseline_time_max_ms 	= 0

opts.each do |opt, arg|
  	case opt
 		when '--trial_amplitude_plot' 		then trial_amplitude_plot	= arg
  		when '--accepted_trial_value_files' 	then accepted_trial_value_files	= arg
 		when '--filelist_filename' 		then filelist_filename 		= arg
		when '--prepended_id_list'   		then prepended_id_list 		= arg
		when '--norm_baseline_type_to_use' 	then norm_baseline_type_to_use 	= arg	
		when '--case_number' 			then case_number 		= arg
		when '--elec_list_type' 		then elec_list_type 		= arg
		when '--use_grid' 			then use_grid 			= arg	
		when '--hp_filter' 			then hp_filter	 		= arg		
		when '--lp_filter' 			then lp_filter	 		= arg		
		when '--pre_stim_time_ms' 		then pre_stim_time_ms 		= arg
		when '--out_min_trials' 		then out_min_trials 		= arg	
		when '--out_max_trials' 		then out_max_trials 		= arg		
		when '--response_window_min_time_ms' 	then response_window_min_time_ms= arg	
		when '--response_window_max_time_ms' 	then response_window_max_time_ms= arg	
		when '--threshold_value' 		then threshold_value 		= arg
		when '--threshold_electrodes' 		then threshold_electrodes	= arg
		when '--threshold_min_time_ms' 		then threshold_min_time_ms 	= arg	
		when '--threshold_max_time_ms' 		then threshold_max_time_ms 	= arg	
		when '--st_baseline_time_min_ms' 	then st_baseline_time_min_ms 	= arg	
		when '--st_baseline_time_max_ms' 	then st_baseline_time_max_ms 	= arg						
  	end
end

t = Time.now
number = t.strftime("%Y%m%d_%H%M%S")

the_working_directory = Dir.getwd
logfile               = "calc_hdf1_st_v60_rb_#{number}.log"

File.delete("#{logfile}") if FileTest.exist?("#{logfile}") 


#os_id = 0
#system("uname -r | cut -f1 -d. > os.var")
#File.open("os.var", "r") do |aFile|
#	aFile.each_line {|line| os_id = line.to_i if (line.strip.length > 0)}
#end 
#File.delete("os.var") if FileTest.exist?("os.var")

if (ARGV.empty? == TRUE) then
	if (filelist_filename.strip.length == 0) then
		puts ""
		print usage
		puts ""
		exit
	end
else
	filelist_filename = "hdf1_files_#{number}.lst"
	File.delete("#{filelist_filename}") if FileTest.exist?("#{filelist_filename}") 
	lstFile = File.new("#{filelist_filename}", "w+")
	count=0
	ARGV.each do |inp_file|
		count = count + 1
		lstFile.puts "#{the_working_directory}/#{inp_file}"
	end
	lstFile.close
end

if FileTest.exist?("#{filelist_filename}") == FALSE then
	puts "missing filelist file #{filelist_filename}: exiting"
	exit
else
	inp_filename_list = IO.readlines("#{filelist_filename}")	
end

#if (os_id.to_i == 2) then
#	if (use_grid.to_i == 1)
#		puts "warning: unable to use grid system from linux os"
#		use_grid = 0
#	end
#end



count_file = 0
inp_filename_list.each do | inp_file_raw |
	count_file    = count_file+1
	hdf_filename = inp_file_raw.strip.chomp
 
      system(" /export/home/nmanz/programs/bourne-shell/calc_hdf1_check_v60.sh #{hdf_filename} #{case_number}")

		mat_files = "mat-file_exists.lst"

       	if FileTest.exist?("#{mat_files}") == TRUE then
				File.delete("#{mat_files}")
			else
       	 #if FileTest.exist?("#{mat_files}") == FALSE then
		puts "Calculating the v6.0 mat file for #{hdf_filename} with e#{elec_list_type} and case #{case_number}"

		replaceFileText("hdf_filename", 						"#{hdf_filename}",     	         	"/export/home/nmanz/programs/matlab/calc_hdf1_st_v6.0.skel", "calc_hdf1_st_v60_#{number}.a")	
		replaceFileText("case_number", 						"#{case_number}", 					 	"calc_hdf1_st_v60_#{number}.a", "calc_hdf1_st_v60_#{number}.b")
		replaceFileText("elec_list_type", 					"#{elec_list_type}", 				 	"calc_hdf1_st_v60_#{number}.b", "calc_hdf1_st_v60_#{number}.a")
		replaceFileText("pre_stim_time_ms", 	        	"#{pre_stim_time_ms}", 	  			 	"calc_hdf1_st_v60_#{number}.a", "calc_hdf1_st_v60_#{number}.b")
		replaceFileText("out_min_trials", 					"#{out_min_trials}", 				 	"calc_hdf1_st_v60_#{number}.b", "calc_hdf1_st_v60_#{number}.a")
		replaceFileText("out_max_trials", 					"#{out_max_trials}", 				 	"calc_hdf1_st_v60_#{number}.a", "calc_hdf1_st_v60_#{number}.b")
		replaceFileText("threshold_value", 					"#{threshold_value}", 	  	 			"calc_hdf1_st_v60_#{number}.b", "calc_hdf1_st_v60_#{number}.a")
		replaceFileText("threshold_min_time_ms", 			"#{threshold_min_time_ms}",   		"calc_hdf1_st_v60_#{number}.a", "calc_hdf1_st_v60_#{number}.b")
		replaceFileText("threshold_max_time_ms", 			"#{threshold_max_time_ms}",   		"calc_hdf1_st_v60_#{number}.b", "calc_hdf1_st_v60_#{number}.a")
		replaceFileText("the_working_directory", 			"#{the_working_directory}",   		"calc_hdf1_st_v60_#{number}.a", "calc_hdf1_st_v60_#{number}.b")
		replaceFileText("norm_baseline_type_to_use", 	"#{norm_baseline_type_to_use}",  	"calc_hdf1_st_v60_#{number}.b", "calc_hdf1_st_v60_#{number}.a")
		replaceFileText("response_window_min_time_ms", 	"#{response_window_min_time_ms}",	"calc_hdf1_st_v60_#{number}.a", "calc_hdf1_st_v60_#{number}.b")
		replaceFileText("response_window_max_time_ms", 	"#{response_window_max_time_ms}",	"calc_hdf1_st_v60_#{number}.b", "calc_hdf1_st_v60_#{number}.a")
		replaceFileText("st_baseline_time_min_ms", 		"#{st_baseline_time_min_ms}", 		"calc_hdf1_st_v60_#{number}.a", "calc_hdf1_st_v60_#{number}.b")
		replaceFileText("st_baseline_time_max_ms", 		"#{st_baseline_time_max_ms}", 		"calc_hdf1_st_v60_#{number}.b", "calc_hdf1_st_v60_#{number}.a")
		replaceFileText("hp_filter", 							"#{hp_filter}", 		 					"calc_hdf1_st_v60_#{number}.a", "calc_hdf1_st_v60_#{number}.b")
		replaceFileText("lp_filter", 							"#{lp_filter}", 		 					"calc_hdf1_st_v60_#{number}.b", "calc_hdf1_st_v60_#{number}.a")
		replaceFileText("threshold_electrodes", 			"#{threshold_electrodes}", 	 		"calc_hdf1_st_v60_#{number}.a", "calc_hdf1_st_v60_#{number}.b")
		replaceFileText("trial_amplitude_plot", 			"#{trial_amplitude_plot}", 		 	"calc_hdf1_st_v60_#{number}.b", "calc_hdf1_st_v60_#{number}.a")
		replaceFileText("accepted_trial_value_files", 	"#{accepted_trial_value_files}", 	"calc_hdf1_st_v60_#{number}.a", "calc_hdf1_st_v60_#{number}.m")
	
		File.delete("calc_hdf1_st_v60_#{number}.a") if FileTest.exist?("calc_hdf1_st_v60_#{number}.a") 
		File.delete("calc_hdf1_st_v60_#{number}.b") if FileTest.exist?("calc_hdf1_st_v60_#{number}.b")

		if (use_grid.to_i == 1) then
			File.copy("#{the_working_directory}/calc_hdf1_st_v60_#{number}.m", "#{the_working_directory}/calc_hdf1_st_v60_#{number}_#{count_file}.m")
			system("grid /opt/bin/matbat10nm -- #{the_working_directory}/calc_hdf1_st_v60_#{number}_#{count_file}.m")
		else
			system("export LD_LIBRARY_PATH=/active_projects/programs/hdf1progs/lib-ibm-linux; /opt/bin/matbat10nm calc_hdf1_st_v60_#{number}.m 2>&1 > #{logfile}")
		end

		File.delete("calc_hdf1_st_v60_#{number}.log") if FileTest.exist?("calc_hdf1_st_v60_#{number}.log")
		File.delete("calc_hdf1_st_v60_#{number}.m")   if FileTest.exist?("calc_hdf1_st_v60_#{number}.m")
	end
end
