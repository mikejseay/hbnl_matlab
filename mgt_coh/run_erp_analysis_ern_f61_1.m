    load('A:\matlab\mgt_coh\subjectid_table_ctl.mat')
    output_dir = 'C:\Users\Mike\Desktop\localdata\active_projects\mike\ern_phase4_ctl';
	data_file_type = 'hdf_binary';
	info_file = [];
	load ern_params_07_15.mat
	param_struct = ern_params;
	option = [true, false, false, false];
	suffix = '32c_20f_18p';
	log_file = 'C:\Users\Mike\Desktop\localdata\active_projects\mike\ern_phase4_ctl\SUNY.phase4_baseline_ern.files.part_1.err';
	if nargin < 8 || isempty(log_file)
		fid = 1;
	else
		fid = fopen(log_file, 'w');
	end
	%filename = textread(file_list, '%s');
    filename = good_list;
	filename = char(filename);
	[n_files, x] = size(filename);
	for m = 1:n_files
		data_file = deblank(filename(m, :));
		if ~exist(data_file, 'file')
			fprintf(fid, 'not found: %s\n', data_file);
			continue
		end
		base_file = basename(data_file);
		output_file = sprintf('%s/%s.%s.mat', output_dir, base_file, suffix);
		if exist(output_file, 'file')
			fprintf(fid, 'skipping %s\n', output_file);
			continue
		else
			fprintf(fid, 'starting %s ',base_file );
            % add a line using touch
            % cmd=sprintf('touch %s', output_file);
            % unix(cmd)
            % or just create an empty file?
		end
		try
			output = erp_analysis(data_file, data_file_type, info_file, ...
				param_struct, option);
		catch err
			fprintf(fid, '%s: ', err.message);
			output = [];
		end
		if ~isempty(output)
			n_trials = sum(output.trial_mat);
			fprintf(fid, '%d ', n_trials);
			fprintf(fid, '\n');
			save(output_file, '-struct', 'output', '-v6');
		else
			fprintf(fid, 'error; not saved\n');
		end
	end
	if fid > 1
		fclose(fid);
	end
