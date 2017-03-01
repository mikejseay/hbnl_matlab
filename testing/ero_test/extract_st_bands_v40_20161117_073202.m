%
% function extract_st_bands_v4.0
%
% modified by Niklas Manz on (2007-06-18)
%
% ln calculation corrected (2008-02-06)
%
% add_info removed (2012-03-14)
% run info added (2012-03-14)
% added the n_chans check at line 245 (2013-03-07)

clear;

path(path,'/export/home/kevjones/progs/hdf1progs/hdf1tools/hdf1matlab/');
path(path,'/export/home/nmanz/programs/matlab/');
path(path,'/export/home/mike/programs/matlab/');

current_time = datestr(now, 'mm-dd-yyyy-HH-MM-SS-FFF');

slog_dir = '/active_projects/ERO_scripts/matlab_slogs/';
slog_filenm = [current_time, '.slog'];
slog_path = fullfile(slog_dir, slog_filenm);
slog_fid = fopen(slog_path, 'w');

try

    do_baseline 				= 0;							% -b (0=no default, 1=yes)
    st_type 						= 1;							% -d (1=total, 2=evoked, 3=induced (tot-evk))
    channel_sort				= 4;					% -e (1=old elec_list default, 2=rows across head, 3=rows within 6 regions, own elec_list)
    electrodes_string 		= '13 16 2 10 18 5 4 19 11';						% -e
    electrodes_array 			= sscanf(electrodes_string, '%f');	% -e
    inp_files 					= '/processed_data/EROprc_lists/4-e1-n10-s9-t100-v800-vp3-tt-64_mats-1479385001913_L14369.lst';				% -f
    calc_type 					= 1;							% -m (1=mean default, 2=max, 3=centroid, 4=maxfreq, 5=maxtime, 6=sum)
    output_text	      	     	= '';						% -o
    out_type 					= 1;							% -p (1=power (lin), 2=amplitude (lin), 3=power (ln), 4=amplitude (ln))
    add_baseline 				= 0;					% -q (0=no, 1=yes default)
    frequencies_min_string 	= '1';				% -v (1= all bands default, 2=single frequencies, 3=lower bands, 4=higher bands, own freq_file)	
    frequencies_max_string 	= '1';				% -v
    freqs_min_array 			= sscanf(frequencies_min_string, '%f');
    freqs_max_array 			= sscanf(frequencies_max_string, '%f');
    n_freq_files 				= min(min([length(freqs_min_array) length(freqs_max_array)]));
    t_min 						= 300;						% -x (300 default)
    t_max 						= 500;						% -y (500 default)

    if     (st_type == 1)
            st_type_name = 'tot';
    elseif (st_type == 2)
            st_type_name = 'evo';
    elseif (st_type == 3)
            st_type_name = 'ind';
    else
            st_type_name = 'tot';
    end

    ln_calc = 0;
    if     (out_type == 1)
            out_type_name = 'pwr';
    elseif (out_type == 2)
            out_type_name    = 'amp';
    elseif (out_type == 3)
            out_type_name = 'pwr-ln';
            out_type = 1;
            ln_calc  = 1;
    elseif (out_type == 4)
            out_type_name = 'amp-ln';
            out_type = 2;
            ln_calc  = 1;
    else
            out_type_name = 'pwr';
    end

    if (length(output_text) == 0)
        out_text = ['v4.0.csv'];
    else
        out_text = ['_v4.0_',output_text,'.csv'];
    end

    [file_names] = textread(inp_files,'%s');
    n_input_files = size(file_names,1);

    %% load first file to get general info
    grand_mean_data = 0;
    load(file_names{1});

    if (exist('h1_gmn_struct'))
        h1_struct       = h1_gmn_struct;
        grand_mean_data = 1;
    end

    %%% NEW %%%

    prc_version = '4';
    exp_name   = h1_struct.st_struct.exp_name;
    case_name  = h1_struct.st_struct.exp_case_type;
    elec_array = strvcat(h1_struct.st_struct.channel_list);
    n_chans    = size(elec_array, 1);

    log_dir = '/active_projects/ERO_scripts/matlab_logs/';
    log_filenm = [prc_version, '-', exp_name, '-', case_name, '-', num2str(n_chans), '-', current_time, '.log'];
    log_path = fullfile(log_dir, log_filenm);
    log_fid = fopen(log_path, 'w');

catch err

    fprintf(slog_fid, '%s: ', err.message);

end


try

    if (grand_mean_data == 1)
        disp 'Warning: program is not intended for grand mean files';
    end

    % ----- inserted to make the following changes for case names of ans experiments:
    if (exp_name == 'ans')
        if (h1_struct.st_struct.exp_case_type == 't1' )
            h1_struct.st_struct.exp_case_type = 'r1';
        end
        if (h1_struct.st_struct.exp_case_type == 't2' )
            h1_struct.st_struct.exp_case_type = 'f2';
        end
        if (h1_struct.st_struct.exp_case_type == 't3' )
            h1_struct.st_struct.exp_case_type = 'f3';
        end
        if (h1_struct.st_struct.exp_case_type == 't4' )
            h1_struct.st_struct.exp_case_type = 'f4';
        end
        if (h1_struct.st_struct.exp_case_type == 't5' )
            h1_struct.st_struct.exp_case_type = 'f5';
        end
        if (h1_struct.st_struct.exp_case_type == 't6' )
            h1_struct.st_struct.exp_case_type = 'f6';
        end
        if (h1_struct.st_struct.exp_case_type == 't7' )
            h1_struct.st_struct.exp_case_type = 'f7';
        end
        if (h1_struct.st_struct.exp_case_type == 't8' )
            h1_struct.st_struct.exp_case_type = 'f8';
        end
        if (h1_struct.st_struct.exp_case_type == 't9' )
            h1_struct.st_struct.exp_case_type = 'f9';
        end
       case_name  = h1_struct.st_struct.exp_case_type;
    end
    % ----- inserted to make the following changes for case names of aod experiments:
    if (exp_name == 'aod')
        if (h1_struct.st_struct.exp_case_type == 't' )
            h1_struct.st_struct.exp_case_type = 'tt';
        end
       case_name  = h1_struct.st_struct.exp_case_type;
    end
    % ----- inserted to make the following changes for case names of cpt experiments:
    if (exp_name == 'cpt')
        if (h1_struct.st_struct.exp_case_type == 't' )
            h1_struct.st_struct.exp_case_type = 'cg';
        end
        if (h1_struct.st_struct.exp_case_type == 'n' )
            h1_struct.st_struct.exp_case_type = 'dn';
        end
        if (h1_struct.st_struct.exp_case_type == 'ng' )
            h1_struct.st_struct.exp_case_type = 'un';
        end
       case_name  = h1_struct.st_struct.exp_case_type;
    end
    % ----- inserted to make the following changes for case names of vp3 experiments:
    if (exp_name == 'vp3')
        if (h1_struct.st_struct.exp_case_type == 't' )
            h1_struct.st_struct.exp_case_type = 'tt';
        end
        if (h1_struct.st_struct.exp_case_type == 'n' )
            h1_struct.st_struct.exp_case_type = 'nv';
        end
       case_name  = h1_struct.st_struct.exp_case_type;
    end


    for k_freq_file=1:n_freq_files

        band1 = freqs_min_array(k_freq_file);
        band2 = freqs_max_array(k_freq_file);
        
         if ( size(findstr(num2str(band1), '.'),1) == 0 )
         	band1 = [num2str(band1),'.0'];
         else
             band1 = num2str(band1);
         end
         if ( size(findstr(num2str(band2), '.'),1) == 0 )
         	band2 = [num2str(band2),'.0'];
         else
             band2 = num2str(band2);
         end

        if (size(case_name,2) == 0)
            out_file = [exp_name,       '_', ...
                band1, '-', ...
                band2, '_', ...
                num2str(t_min), '-', ...
                num2str(t_max), ...
                '_b', num2str(do_baseline), ...
                '_m', num2str(calc_type), '_', ...
                st_type_name, out_text];
            out_base_file = [exp_name, '_', band1, '-', band2, ...
                '_base', '_m', num2str(calc_type), '_', st_type_name,'-', out_type_name, out_text];
        else
            out_file = [exp_name, '-', case_name, '_', band1, '-', band2, '_', num2str(t_min), '-', num2str(t_max), ...
                '_b', num2str(do_baseline), '_m', num2str(calc_type), '_', st_type_name,'-', out_type_name, out_text];
            out_base_file = [exp_name, '-', case_name, '_', band1, '-', band2, ...
                '_base', '_m', num2str(calc_type), '_', st_type_name,'-', out_type_name, out_text];
        end

        if (channel_sort == 2) %sorting along the rows across the head
            if (n_chans == 19)
                new_sort = [15,17, 12,13,16,2,3, 8,10,18,5,7, 6,4,19,11,9, 1,14];
            elseif (n_chans == 20)
                n_chans = 19;
                new_sort = [15,17, 12,13,16,2,3, 8,10,18,5,7, 6,4,19,11,9, 1,14];
            elseif (n_chans == 32)
                n_chans  = 31;
                new_sort = [1,2, 5,6, 3,9,7,8,4, 11,13,12,10, 15,17,16,18,14, 19,21,22,20, 27,23,25,24,26, 29,28, 31,30];
            elseif (n_chans == 62)
                n_chans  = 61;
                new_sort = [1,39,2, 33,5,48,6,34, 3,35,9,45,7,44,8,36,4, 37,11,41,13,57,12,40,10,38, 15,43,17,53,16,54,18,42,14, 47,19,49,21,62,22,50,20,46, 27,51,23,61,25,60,24,52,26, 55,29,58,28,56, 31,59,30];
            end

        end
        if (channel_sort == 3) %sorting along the rows in 6 regions
            if (n_chans == 19)
                new_sort = [15,17, 12,13,16,2,3, 10,18,5, 4,19,11, 1,14,  8,6, 7,9];
            elseif (n_chans == 20)
                n_chans = 19;
                new_sort = [15,17, 12,13,16,2,3, 8,10,18,5,7, 6,4,19,11,9, 1,14];
            elseif (n_chans == 32)
                n_chans  = 31;
                new_sort = [1,2,5,6, 3,9,7,8,4, 13,12, 17,16,18, 21,22,23,25,24, 29,28,31,30, 11,15,19,27, 10,14,20,26];
            elseif (n_chans == 62)
                n_chans  = 61;
                new_sort = [1,39,2, 33,5,48,6,34, 3,35,9,45,7,44,8,36,4, 11,41,13,57,12,40,10, 17,53,16,54,18, 49,21,62,22,50, 23,61,25,60,24, 55,29,58,28,56, 31,59,30, 37,15,43,47,19,27,51, 38,42,14,20,46,52,26];
            end

        end

    end

    for i_file = 1:n_input_files

        try

            filenm = file_names{i_file};

            %%% NEW %%%
            try
                load(filenm); % replaced by block
            catch
                fprintf(log_fid, '%s was missing\n', filenm)
                continue
            end
            %%% NEW %%%

            grand_mean_data = 0;
            if (exist('h1_gmn_struct'))
                h1_struct = h1_gmn_struct;
                grand_mean_data = 1;
            end

            if (grand_mean_data == 1)
                file_id = '99999999';
                file_session = 'z';
                file_run     = '9';

                disp 'Warning: program is not intended for grand mean files';
                continue
            else
                file_id      = h1_struct.subject_struct.subject_id{1};
               
                if ( file_id(2) == '5' )
                    file_id(2) = '0';
                    h1_struct.subject_struct.run_number{1} = '2';
                end    
                file_session = h1_struct.subject_struct.session_code{1};
                file_run     = h1_struct.subject_struct.run_number{1};
                if     ( file_id(1) == 'a' ) & ( file_id(1:5) ~= 'a0000')
                         file_id(1) = '1';
                elseif ( file_id(1) == 'b' )
                         file_id(1) = '2';
                elseif ( file_id(1) == 'c' ) & ( file_id(1:5) ~= 'c0000')
                         file_id(1) = '3';
                elseif ( file_id(1) == 'd' )
                         file_id(1) = '4';
                elseif ( file_id(1) == 'e' )
                         file_id(1) = '5';
                elseif ( file_id(1) == 'f' )
                         file_id(1) = '6';
                end    
            end

        %    %%%%%%%%%%%%%%%%%%%%%%%%
        %    % cutting h1_struct.st_struct.n_chans down if size(h1_struct.data_struct.hdf1_avg_st_data) is smaller
        %
        %    if ( size(h1_struct.data_struct.hdf1_avg_st_data,1) < n_chans )
        %         n_chans = size(h1_struct.data_struct.hdf1_avg_st_data,1);
        %    end
        %    %%%%%%%%%%%%%%%%%%%%%%%%

            %%% NEW CODE FOR SAVING 3D ARRAYS OF S_VALUES %%%

            % clip times before entering into s-transform so that frequencies
            % are uniformly chosen

            % hardcode_start_ms = -187.5;
            % hardcode_end_ms = 1000;
            
            % if abs(hardcode_start_ms) > h1_struct.experiment_struct.pre_stim_time_ms
            %     fprintf(log_fid, '%s did not have enough pre-stimulus time\n', filenm)
            %     continue
            % end

            % if hardcode_end_ms > h1_struct.experiment_struct.post_stim_time_ms
            %     fprintf(log_fid, '%s did not have enough post-stimulus time\n', filenm)
            %     continue
            % end
            
            % rate = h1_struct.experiment_struct.rate;
            % time_vec_tmp = (-h1_struct.experiment_struct.pre_stim_time_ms:(1000/rate): ...
            %     h1_struct.experiment_struct.post_stim_time_ms)';
            % [junk, start_idx] = min(abs(time_vec_tmp - hardcode_start_ms));
            % [junk, end_idx] = min(abs(time_vec_tmp - hardcode_end_ms));

            % h1_struct.data_struct.hdf1_avg_st_data = h1_struct.data_struct.hdf1_avg_st_data(:, start_idx:end_idx);
            % h1_struct.data_struct.hdf1_ind_st_data = h1_struct.data_struct.hdf1_ind_st_data(:, start_idx:end_idx);

            [S_values, Sbase_values, S_all, t, f] = calc_chanlist_stockwell_values_ms(h1_struct, ...
                freqs_min_array, freqs_max_array, t_min, t_max, calc_type, ...
                h1_struct.st_struct.channel_list, st_type, do_baseline, ...
                h1_struct.st_struct.baseline_type, h1_struct.st_struct.st_baseline_time_min_ms, ...
                h1_struct.st_struct.st_baseline_time_max_ms, out_type);

            f_ds_vec = (f>=0.5 & f<=50)';
            f_vec = f';
            f_vec_ds = f(f_ds_vec)';
            
            rate = h1_struct.experiment_struct.rate;
            time_ds_factor = 2;
            time_vec = (-h1_struct.experiment_struct.pre_stim_time_ms:(1000/rate): ...
            h1_struct.experiment_struct.post_stim_time_ms)';
            time_cut_vec = (time_vec>=-400 & time_vec<=1000);
            time_downsamp_vec = false(size(time_vec));
            time_downsamp_vec(1:time_ds_factor:end) = true;
            time_ds_vec = time_cut_vec & time_downsamp_vec;
            time_vec_ds = time_vec(time_ds_vec);

            % time_ds_factor = 2;
            % time_vec = (-hardcode_start_ms:(1000/rate):hardcode_end_ms)';
            % time_ds_vec = false(size(time_vec));
            % time_ds_vec(1:time_ds_factor:end) = true;
            % time_vec_ds = time_vec(time_ds_vec);

            n_chans_present = size(S_all, 3);
            
            n_trials = h1_struct.st_struct.num_trials;
            st_baseline_time_min_ms = h1_struct.st_struct.st_baseline_time_min_ms;
            st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
            hp_filter = h1_struct.st_struct.hp_filter;
            lp_filter = h1_struct.st_struct.lp_filter;
            
            opt = v2struct(add_baseline, calc_type, case_name, channel_sort, do_baseline, ...
                elec_array, exp_name, file_id, filenm, file_run, ...
                file_session, ln_calc, n_chans, n_chans_present, ...
                out_type, out_type_name, ...
                st_type, st_type_name, ...
                f_vec, f_vec_ds, rate, time_ds_factor, time_vec, time_vec_ds, ...
                n_trials, st_baseline_time_min_ms, st_baseline_time_max_ms, hp_filter, lp_filter);
            
            S_all_ds = S_all(f_ds_vec, time_ds_vec, :);
            clear S_all
            data = single(S_all_ds);
            clear S_all_ds
            
            folderparts = strread(filenm,'%s','delimiter','/');
            nameparts = strread(folderparts{end},'%s','delimiter','_');
            session_part = nameparts{3};
            id = nameparts{4};
            session = session_part(1);
            param_str = folderparts{end-2};
            
            parent_dir = '/processed_data/ero-mats-v40';
            outname = [strjoin({id, session, exp_name, case_name, st_type_name}, '_'),'.mat'];
            outdir = fullfile(parent_dir, param_str, num2str(n_chans_present), exp_name);
            outpath = fullfile(outdir, outname);
            
            if ~exist(outdir, 'dir') %make if DNE yet
                mkdir(outdir);
            end
            
            if ~exist(outpath, 'file')
                save(outpath, 'data', 'opt', '-v7.3');
            end
            %%% NEW CODE FOR SAVING 3D ARRAYS OF S_VALUES %%%

        catch err
            fprintf(log_fid, '%s: ', err.message);
            continue
        end

    end

catch err

    fprintf(log_fid, '%s: ', err.message);

end

fclose(slog_fid);
fclose(log_fid);
