function [h1_struct] = create_cnthdf1_baseline_epoch_data_v6(h1_struct, new_txt_file_name, new_hi_pass, new_lo_pass, new_threshold_electrodes, the_case, trial_plot)

%
% adding new_txt_file_name, new_hi_pass, and new_lo_pass, new_threshold_electrodes (Niklas, 2008-07-22)
%
% correcting a problem with the case_num (if larger than size(case_num,2)) in line 251

fid = fopen(new_txt_file_name, 'a');

%scale the cnt data

if (h1_struct.file_struct.data_type ~= 0)
    return
end

% set number of electrodes for threshold checking

fprintf(fid, 'Number of electrodes           : %d \n', h1_struct.experiment_struct.n_chans);
if (h1_struct.experiment_struct.n_chans == 21) 
   fprintf(fid, 'Number of head electrodes      : 19\n');
elseif (h1_struct.experiment_struct.n_chans == 32) 
   fprintf(fid, 'Number of head electrodes      : 31\n');
elseif (h1_struct.experiment_struct.n_chans == 64) 
   fprintf(fid, 'Number of head electrodes      : 61\n');
end

all_electr = 0;
if (size(new_threshold_electrodes,2) == 1 )
    if ( new_threshold_electrodes(1) == 61 )
        all_electr = 1;
    end
end
if (h1_struct.experiment_struct.n_chans == 21) && (size(new_threshold_electrodes,2) == 0)
    new_threshold_electrodes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19];
    n_threshold_chans = 19;
    fprintf(fid, 'Number of threshold electrodes : %2.0f \n', n_threshold_chans);
    fprintf(fid, 'Which threshold electrodes     : all 19 head electrodes \n');
    for elec=1:size(new_threshold_electrodes,2)
        if ( new_threshold_electrodes(elec) < 3 )
            new_threshold_electrode_number(elec) = new_threshold_electrodes(elec);
        elseif   ( new_threshold_electrodes(elec) >= 3 ) && ( new_threshold_electrodes(elec) < 13 )
            new_threshold_electrode_number(elec) = new_threshold_electrodes(elec)+1;
        elseif  ( new_threshold_electrodes(elec) >= 13 )
            new_threshold_electrode_number(elec) = new_threshold_electrodes(elec)+2;
        end
    end
    elec_thresh_names = h1_struct.run_struct.channel_label(new_threshold_electrode_number);
elseif (h1_struct.experiment_struct.n_chans == 21) && (size(new_threshold_electrodes,2) > 1)
    n_threshold_chans = size(new_threshold_electrodes,2);
    fprintf(fid, 'Number of threshold electrodes : %2.0f \n', n_threshold_chans);
    fprintf(fid, 'Which threshold electrodes     : ');
    for elec=1:size(new_threshold_electrodes,2)
        if ( new_threshold_electrodes(elec) < 3 )
            new_threshold_electrode_number(elec) = new_threshold_electrodes(elec);
        elseif   ( new_threshold_electrodes(elec) >= 3 ) && ( new_threshold_electrodes(elec) < 13 )
            new_threshold_electrode_number(elec) = new_threshold_electrodes(elec)+1;
        elseif  ( new_threshold_electrodes(elec) >= 13 )
            new_threshold_electrode_number(elec) = new_threshold_electrodes(elec)+2;
        end
    end
    elec_thresh_names = h1_struct.run_struct.channel_label(new_threshold_electrode_number);
    for i=1:size(new_threshold_electrodes,2)-1
        fprintf(fid, '%s(%s), ' ,num2str(new_threshold_electrodes(i)),elec_thresh_names{i});
    end
    fprintf(fid, '%s(%s)\n' ,num2str(new_threshold_electrodes(size(new_threshold_electrodes,2))),elec_thresh_names{size(new_threshold_electrodes,2)});
elseif (h1_struct.experiment_struct.n_chans == 32) && (size(new_threshold_electrodes,2) == 0)
    new_threshold_electrodes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31];
    n_threshold_chans = 31;
    fprintf(fid, 'Number of threshold electrodes : %2.0f \n', n_threshold_chans);
    fprintf(fid, 'Which threshold electrodes     : all 31 head electrodes\n');
elseif (h1_struct.experiment_struct.n_chans > 32) && (size(new_threshold_electrodes,2) == 0)
    new_threshold_electrodes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31];
    n_threshold_chans = 31;
    fprintf(fid, 'Number of threshold electrodes : %2.0f \n', n_threshold_chans);
    fprintf(fid, 'Which threshold electrodes     : first 31 head electrodes\n');
elseif (h1_struct.experiment_struct.n_chans > 32) && ( all_electr == 1 )
    new_threshold_electrodes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62];
    n_threshold_chans = 61;
    fprintf(fid, 'Number of threshold electrodes : %2.0f \n', n_threshold_chans);
    fprintf(fid, 'Which threshold electrodes     : all 61 head electrodes\n');
elseif (h1_struct.experiment_struct.n_chans > 32) && (size(new_threshold_electrodes,2) == 1)
    new_threshold_electrodes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31];
    n_threshold_chans = 31;
    fprintf(fid, 'Number of threshold electrodes : %2.0f \n', n_threshold_chans);
    fprintf(fid, 'Which threshold electrodes     : first 31 head electrodes\n');
elseif (h1_struct.experiment_struct.n_chans >= 32) && (size(new_threshold_electrodes,2) > 1)
    n_threshold_chans = size(new_threshold_electrodes,2);
    fprintf(fid, 'Number of threshold electrodes : %2.0f \n', n_threshold_chans);
    fprintf(fid, 'Which threshold electrodes     : ');
    elec_thresh_names = h1_struct.run_struct.channel_label(new_threshold_electrodes);
    for i=1:size(new_threshold_electrodes,2)-1
        fprintf(fid, '%s(%s), ' ,num2str(new_threshold_electrodes(i)),elec_thresh_names{i});
    end
    fprintf(fid, '%s(%s)\n' ,num2str(new_threshold_electrodes(size(new_threshold_electrodes,2))),elec_thresh_names{size(new_threshold_electrodes,2)});
end
if ( h1_struct.experiment_struct.n_chans == 21 )
 new_threshold_electrodes    = new_threshold_electrode_number;
end


fprintf(fid, '\n');

rate                 = double(h1_struct.experiment_struct.rate);
trial_time_length_ms = h1_struct.experiment_struct.pre_stim_time_ms + h1_struct.experiment_struct.post_stim_time_ms;
num_trial_samples    = floor(trial_time_length_ms / ((1.0/rate) * 1000.));
pre_stim_samples     = floor(h1_struct.experiment_struct.pre_stim_time_ms / ((1.0/rate) * 1000.));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate trial baseline parameters

baseline_start_time_ms = h1_struct.experiment_struct.baseline_time_min_ms;
baseline_stop_time_ms  = h1_struct.experiment_struct.baseline_time_max_ms;
baseline_start_sample  = floor(baseline_start_time_ms / ((1.0/rate) * 1000.)) + pre_stim_samples + 1;
baseline_stop_sample   = floor(baseline_stop_time_ms  / ((1.0/rate) * 1000.)) + pre_stim_samples + 1;
num_baseline_samps     = baseline_stop_sample - baseline_start_sample + 1;

if (baseline_start_sample <= 0)
    baseline_start_sample = 1;
end
if (baseline_stop_sample > num_trial_samples)
    baseline_stop_sample = num_trial_samples;
end

% calculate trial threshold parameters

threshold_start_time_ms = h1_struct.experiment_struct.threshold_time_min_ms;
threshold_stop_time_ms  = h1_struct.experiment_struct.threshold_time_max_ms;
threshold_start_sample  = floor(threshold_start_time_ms / ((1.0/rate) * 1000.)) + pre_stim_samples + 1;
threshold_stop_sample   = floor(threshold_stop_time_ms  / ((1.0/rate) * 1000.)) + pre_stim_samples + 1;
num_threshold_samps     = threshold_stop_sample - threshold_start_sample + 1 ;

if (threshold_start_sample <= 0)
    threshold_start_sample = 1;
end
if (threshold_stop_sample > num_trial_samples)
    threshold_stop_sample = num_trial_samples;
end

% Renaming some case names for the cpt and vp3 experiment.
if ( h1_struct.experiment_struct.exp_name == 'cpt')
    for i = 1:length(h1_struct.case_struct.case_type)
        if (h1_struct.case_struct.case_type{i} == 'T' )
            h1_struct.case_struct.case_type{i}     = 'CG';
        end
        if (h1_struct.case_struct.case_type{i} == 'NG' )
            h1_struct.case_struct.case_type{i}     = 'UN';
        end
        if (h1_struct.case_struct.case_type{i} == 'N' )
            h1_struct.case_struct.case_type{i}     = 'DN';
        end
    end
end
if ( h1_struct.experiment_struct.exp_name == 'vp3')
    for i = 1:length(h1_struct.case_struct.case_type)
        if (h1_struct.case_struct.case_type{i} == 'T' )
            h1_struct.case_struct.case_type{i}     = 'TT';
        end
        if (h1_struct.case_struct.case_type{i} == 'N' )
            h1_struct.case_struct.case_type{i}     = 'NV';
        end
    end
end

% determine filter raw data parameters

if (h1_struct.transforms_struct.lo_pass_filter == 0.0)
    h1_struct.transforms_struct.lo_pass_filter = 55.0;
end

if (h1_struct.transforms_struct.hi_pass_filter == 0.0)
    h1_struct.transforms_struct.hi_pass_filter = 0.05;
end

if (new_lo_pass > 0.0)
    h1_struct.transforms_struct.lo_pass_filter = new_lo_pass;
end

if (new_hi_pass > 0.0)
    h1_struct.transforms_struct.hi_pass_filter = new_hi_pass;
end

if (h1_struct.transforms_struct.lo_pass_filter > 0.0)
    Wn_lo     = h1_struct.transforms_struct.lo_pass_filter / (rate/2);
    %[blo,alo] = cheby1(5,0.5,Wn_lo,'low');
    [blo,alo] = cheby1(5,0.043,Wn_lo,'low');
end

if (h1_struct.transforms_struct.hi_pass_filter > 0.0)
    Wn_hi     = h1_struct.transforms_struct.hi_pass_filter / (rate/2);
    %[bhi,ahi] = cheby1(3,0.5,Wn_hi,'high');
    [bhi,ahi] = cheby1(5,0.043,Wn_hi,'high');
end


%%%%%%%%%%%%  part to remove amplitude values above 10 mV and set them to the previous sample value (Niklas, 2013-05-08)

if ( h1_struct.experiment_struct.uv_per_unit(1) < 0.1 )
    trial = h1_struct.experiment_struct.n_trials-2;
    trial_stimulus_time_ms = h1_struct.trial_struct.time_offset(trial) * 1000.0;
    trial_stimulus_sample  = floor(trial_stimulus_time_ms / ((1.0/rate) * 1000.));
    trial_start_sample     = trial_stimulus_sample - pre_stim_samples + 1;

    position_start = trial_start_sample;
    position_end = size(h1_struct.data_struct.hdf1_cnt_data,2);
    for channel = 1:32
        position = position_start;
        while position < position_end
            samp_value = abs(h1_struct.data_struct.hdf1_cnt_data(channel,position));
            if (samp_value > 10000 ) && (abs(h1_struct.data_struct.hdf1_cnt_data(channel,position+1)) > 10000 )
                for i = 1:h1_struct.experiment_struct.n_chans
                    h1_struct.data_struct.hdf1_cnt_data(i,position-1:size(h1_struct.data_struct.hdf1_cnt_data,2))= h1_struct.data_struct.hdf1_cnt_data(i,position-2);
                end
                position = position_end;
                fprintf(fid, '\n');
                fprintf(fid, 'The experiment has an amplitude value of %7.0f muV at sample point %7.0f for channel %s. \nAll following amplitudes have been set to the previous value.\n ', samp_value, position, num2str(channel));
           end
            position = position + 1;
        end
    end
    fprintf(fid, '\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avg_values = mean (h1_struct.data_struct.hdf1_cnt_data,2);

for channel = 1:h1_struct.experiment_struct.n_chans
    
    % determine average correction value
    trace = (h1_struct.data_struct.hdf1_cnt_data(channel,:) - avg_values(channel));

    % lo pass filter trace data
    if (h1_struct.transforms_struct.lo_pass_filter > 0.0)
        trace = filtfilt(blo,alo,trace);
    end

    % hi pass filter trace data
    if (h1_struct.transforms_struct.hi_pass_filter > 0.0)
        trace = filtfilt(bhi,ahi,trace);
    end

    h1_struct.data_struct.hdf1_cnt_data(channel,:) = trace;

end

h1_struct.data_struct.hdf1_epoch_data = zeros(h1_struct.experiment_struct.n_trials,h1_struct.experiment_struct.n_chans,num_trial_samples);
if (trial_plot == 1)
    h1_struct.data_struct.hdf1_orig_data = zeros(h1_struct.experiment_struct.n_trials,h1_struct.experiment_struct.n_chans,num_trial_samples);
end

n = 0;
y = 0;

for trial = 1:h1_struct.experiment_struct.n_trials-1
    trial_stimulus_time_ms = h1_struct.trial_struct.time_offset(trial) * 1000.0;
    trial_stimulus_sample  = floor(trial_stimulus_time_ms / ((1.0/rate) * 1000.));
    trial_start_sample     = trial_stimulus_sample - pre_stim_samples + 1;
    trial_stop_sample      = trial_start_sample + num_trial_samples;
    %%%%%%%%%% added on 2014-03-04 because of "??? Index exceeds matrix
    %%%%%%%%%% dimensions." trial_stop_sample was larger than size(h1_struct.data_struct.hdf1_cnt_data,2)
    if (size(h1_struct.data_struct.hdf1_cnt_data,2) < trial_stop_sample)
        trial_stop_sample     = size(h1_struct.data_struct.hdf1_cnt_data,2);
    end
    if (size(h1_struct.data_struct.hdf1_cnt_data,2) < threshold_stop_sample)
        threshold_stop_sample = size(h1_struct.data_struct.hdf1_cnt_data,2);
    end
   %%%%%%%%
    for channel = 1:h1_struct.experiment_struct.n_chans
        trial_trace = squeeze(h1_struct.data_struct.hdf1_cnt_data(channel,trial_start_sample:trial_stop_sample));

        % determine baseline value
        if (baseline_start_sample == baseline_stop_sample)
            baseline_value = 0;
        else
            baseline_value = mean(trial_trace(baseline_start_sample:baseline_stop_sample));
        end
        avg_trial_value = mean(trial_trace(:));

        % determine channel max value
        
        electrode_max_value(channel) = max(abs(trial_trace(threshold_start_sample:threshold_stop_sample) - baseline_value));

        % store epoch data
        for i = 1:length(trial_trace)
            if (trial_plot == 1)
                h1_struct.data_struct.hdf1_orig_data(trial,channel,i)  = trial_trace(i);
            end
            h1_struct.data_struct.hdf1_epoch_data(trial,channel,i) = trial_trace(i) - baseline_value;
        end
    end

    % determine threshold artifacts
    if (size(new_threshold_electrodes,2) == 0)
        max_data_value = max(electrode_max_value(1:n_threshold_chans));
        max_data_num   = find(electrode_max_value(1:n_threshold_chans)==max(max(electrode_max_value(1:n_threshold_chans))));
    else
        max_data_value = max(electrode_max_value(new_threshold_electrodes));
        max_data_place = new_threshold_electrodes(find(electrode_max_value(new_threshold_electrodes)==max(max(electrode_max_value(new_threshold_electrodes)))));
        if ( h1_struct.experiment_struct.n_chans == 21 )
            if ( max_data_place < 3 )
                max_data_num = max_data_place;
            elseif   ( max_data_place >= 3 ) && ( max_data_place < 13 )
                max_data_num = max_data_place-1;
            elseif  ( max_data_place >= 13 )
                max_data_num = max_data_place-2;
            end
        else
            max_data_num   = max_data_place;
        end
    end
    elec_name = h1_struct.run_struct.channel_label{max_data_place};
    h1_struct.trial_struct.max_thresh_win_value(trial) = max_data_value;
    h1_struct.trial_struct.threshold_value             = h1_struct.experiment_struct.threshold_value;

    case_position = find(h1_struct.case_struct.case_num==h1_struct.trial_struct.case_num(trial));
    if (case_position<=size(h1_struct.case_struct.case_num,2))
        case_name     = h1_struct.case_struct.case_type{case_position};
    else
        case_name = 'Unknown';
    end

    if (h1_struct.trial_struct.max_thresh_win_value(trial) > h1_struct.experiment_struct.threshold_value)
        h1_struct.trial_struct.artf_present(trial) = 1;
        n = n+1;
        fprintf(fid, 'For trial %3.0f (case %s, %s): above threshold value (%3.0f) at electrode %s(%s)\n', trial, num2str(h1_struct.trial_struct.case_num(trial)), case_name, h1_struct.trial_struct.max_thresh_win_value(trial), num2str(max_data_num),elec_name);
    else
        h1_struct.trial_struct.artf_present(trial) = 0;
        y = y+1;
    end
end
fprintf(fid, '\n');
fprintf(fid, 'Number of trials above threshold: %3.0f \n', n);
fprintf(fid, 'Number of trials below threshold: %3.0f \n', y);
fprintf(fid, '\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addition to plot the trial values of the center 9 electrodes before and after epoch baselining

if (trial_plot == 1)

    if (h1_struct.experiment_struct.n_chans > 20)

        titlestring = replace_underscore(h1_struct.file_struct.file_name(1:17));
        for j = 1:9
            if max(h1_struct.experiment_struct.n_chans) > 21
                if j==1;     k=9;
                elseif j==2; k=7;
                elseif j==3; k=8;
                elseif j==4; k=17;
                elseif j==5; k=16;
                elseif j==6; k=18;
                elseif j==7; k=23;
                elseif j==8; k=25;
                elseif j==9; k=24;
                end
            else
                if j==1;     k=13;
                elseif j==2; k=16;
                elseif j==3; k=2;
                elseif j==4; k=10;
                elseif j==5; k=18;
                elseif j==6; k=5;
                elseif j==7; k=4;
                elseif j==8; k=19;
                elseif j==9; k=11;
                end
            end
            for i = 1:length(h1_struct.trial_struct.trial_num)
                 % baselined values
                maxi=max(h1_struct.data_struct.hdf1_epoch_data(i,k,:));
                mini=min(h1_struct.data_struct.hdf1_epoch_data(i,k,:));
                if abs(mini)>maxi
                    max_total_trial_epoch(j,i)=mini;
                else
                    max_total_trial_epoch(j,i)=maxi;
                end
               % original values
                maxi=max(h1_struct.data_struct.hdf1_orig_data(i,k,:));
                mini=min(h1_struct.data_struct.hdf1_orig_data(i,k,:));
                if abs(mini)>maxi
                    max_total_orig(j,i)=mini;
                else
                    max_total_orig(j,i)=maxi;
                end
            end
        end
        max_orig  = max(max(max_total_orig(:,:)));
        max_epoch = max(max(max_total_trial_epoch(:,:)));

        figure;
        subplot(3,1,1);
        plot(h1_struct.data_struct.hdf1_cnt_data(1,:));
        for j = 2:9
            hold on, plot(h1_struct.data_struct.hdf1_cnt_data(j,:));
        end
        plot([0,h1_struct.data_struct.hdf1_data_dims(2)],[100,100],'r');
        plot([0,h1_struct.data_struct.hdf1_data_dims(2)],[0,0],'r');
        xlim([0, h1_struct.data_struct.hdf1_data_dims(2)]);ylim([-500, 500]);
        xlabel('Sample number'), ylabel('Amplitude (\muV)'), title(['Sample values of the center 9 electrodes for ', titlestring]);

        subplot(3,1,2);
        plot(max_total_orig(1,:));
        for j = 2:9
            hold on, plot(max_total_orig(j,:));
        end
        plot([0,h1_struct.data_struct.hdf1_data_dims(2)],[100,100],'r');
        plot([0,h1_struct.data_struct.hdf1_data_dims(2)],[0,0],'r');
        xlim([0, length(h1_struct.trial_struct.trial_num)]);ylim([-max_orig, max_orig]);
        xlabel('Trial number'), ylabel('Amplitude (\muV)'), title(['Initial max trial values of the center 9 electrodes for ', titlestring]);

        subplot(3,1,3);
        plot(max_total_trial_epoch(1,:));
        for j = 2:9
            hold on, plot(max_total_trial_epoch(j,:));
        end
        plot([0,h1_struct.data_struct.hdf1_data_dims(2)],[100,100],'r');
        plot([0,h1_struct.data_struct.hdf1_data_dims(2)],[0,0],'r');
        xlim([0, length(h1_struct.trial_struct.trial_num)]);ylim([-max_epoch, max_epoch]);
        xlabel('Trial number'), ylabel('Amplitude (\muV)'), title(['Baselined max trial values of the center 9 electrodes for ', titlestring]);

        filestring = h1_struct.file_struct.file_name(1:17);
        plot_file  = [h1_struct.file_struct.file_name(1:17) '_amplitude-info.ps'];
        set(gcf, 'PaperPosition', [0.25, 0.25, 8, 10.5]);
        print(gcf, '-dpsc2', plot_file);
    end

end