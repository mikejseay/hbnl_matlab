function [S, inv_S, t, f] = calc_chan_stockwell_power(h1_struct, chan_name, st_type, apply_baseline, baseline_type, st_baseline_time_min_ms, st_baseline_time_max_ms, out_type)
%
% function [S, inv_S, t, f] = calc_chan_stockwell_power(h1_struct, chan_name, [st_type], [apply_baseline], [baseline_type], [st_baseline_time_min_ms], [st_baseline_time_max_ms], [out_type])
%
if (nargin < 8)
   out_type = 1;
end 
if (nargin < 7)
   out_type = 1;
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
end 
if (nargin < 6)
   out_type = 1; 
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
   st_baseline_time_min_ms = h1_struct.st_struct.st_baseline_time_min_ms;
end 
if (nargin < 5)
   out_type = 1; 
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
   st_baseline_time_min_ms = h1_struct.st_struct.st_baseline_time_min_ms;
   baseline_type = h1_struct.st_struct.baseline_type;
end 
if (nargin < 4)
   out_type = 1; 
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
   st_baseline_time_min_ms = h1_struct.st_struct.st_baseline_time_min_ms;
   baseline_type = h1_struct.st_struct.baseline_type;
   apply_baseline = 0;
end 
if (nargin < 3)
   out_type = 1; 
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
   st_baseline_time_min_ms = h1_struct.st_struct.st_baseline_time_min_ms;
   baseline_type = h1_struct.st_struct.baseline_type;
   apply_baseline = 0;
   st_type = 1;
end 

rate = h1_struct.experiment_struct.rate;
pre_stim_time_ms = h1_struct.experiment_struct.pre_stim_time_ms;

the_case = h1_struct.st_struct.the_case;
case_idx = get_case_idx(h1_struct.case_struct.case_num, the_case);
channel_idx = get_channel_index(h1_struct.st_struct.channel_list, chan_name);

try 
    test_exist = h1_struct.data_struct.hdf1_avg_st_data;
catch
   
   channel_indices = get_channel_indices(h1_struct.run_struct.channel_label, h1_struct.st_struct.channel_list);
 
   if (length(size(h1_struct.data_struct.hdf1_avg_data)) == 3)
       h1_struct.data_struct.hdf1_avg_st_data = squeeze(h1_struct.data_struct.hdf1_avg_data(case_idx,channel_indices,:));
   else
       h1_struct.data_struct.hdf1_avg_st_data = squeeze(h1_struct.data_struct.hdf1_avg_data(channel_indices,:));
   end
   
end 

x_avg = squeeze(h1_struct.data_struct.hdf1_avg_st_data(channel_idx,:));
x_ind = squeeze(h1_struct.data_struct.hdf1_ind_st_data(channel_idx,:));
x_base_ind = squeeze(h1_struct.data_struct.hdf1_ind_st_base_data(channel_idx,:));

[st_x_avg,t,f] = stockwell_transform(x_avg, 0, length(x_avg)/2, 1/rate);
[st_x_ind,t,f] = stockwell_transform(x_ind, 0, length(x_ind)/2, 1/rate);
% [st_x_base_ind,t,f] = stockwell_transform(x_base_ind, 0, length(x_base_ind)/2, 1/rate);


if (apply_baseline == 1)
    st_x_avg_2_complex = apply_tfr_complex_baseline(st_x_avg, st_baseline_time_min_ms, st_baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type);
    abs_base_S = abs(st_x_avg_2_complex);
    angle_base_S = angle(st_x_avg);
    st_x_avg_2_complex = amplitude_phase_to_complex(abs_base_S, angle_base_S);
    
    st_x_ind_2_complex = apply_tfr_complex_baseline(st_x_base_ind, st_baseline_time_min_ms, st_baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type);    
    abs_base_S = abs(st_x_ind_2_complex);
    angle_base_S = angle(st_x_base_ind);
    st_x_ind_2_complex = amplitude_phase_to_complex(abs_base_S, angle_base_S);
    
    log_transform = 0;
    st_x_avg_2_abs = apply_tfr_abs_baseline(st_x_avg, st_baseline_time_min_ms, st_baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type, log_transform);
    st_x_ind_2_abs = apply_tfr_abs_baseline(st_x_base_ind, st_baseline_time_min_ms, st_baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type, log_transform);    

    if (out_type >= 3)
        log_transform = 1;
        st_x_avg_2_log = apply_tfr_abs_baseline(st_x_avg, st_baseline_time_min_ms, st_baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type, log_transform);
        st_x_ind_2_log = apply_tfr_abs_baseline(st_x_base_ind, st_baseline_time_min_ms, st_baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type, log_transform);    
    end
else 
    st_x_avg_2_complex = st_x_avg;
    st_x_ind_2_complex = st_x_ind;
    
    st_x_avg_2_abs = abs(st_x_avg);
    st_x_ind_2_abs = abs(st_x_ind);
    
    if (out_type >= 3)
        st_x_avg_2_log = log(abs(st_x_avg));
        st_x_ind_2_log = log(abs(st_x_ind));
    end
end    

min_f = min([size(st_x_avg_2_complex,1) size(st_x_ind_2_complex,1)]);
min_t = min([size(st_x_avg_2_complex,2) size(st_x_ind_2_complex,2)]);
t = t(1:min_t);
f = f(1:min_f);

if (st_type == 1) % ind
    S_complex = st_x_ind_2_complex(1:min_f,1:min_t);
    S_abs = st_x_ind_2_abs(1:min_f,1:min_t);
    if (out_type >= 3) S_log = st_x_ind_2_log(1:min_f,1:min_t); end;
elseif (st_type == 2) % avg
    S_complex = st_x_avg_2_complex(1:min_f,1:min_t);
    S_abs = st_x_avg_2_abs(1:min_f,1:min_t);
    if (out_type >= 3) S_log = st_x_avg_2_log(1:min_f,1:min_t); end;
elseif (st_type == 3) % ind - avg
    S_complex = st_x_ind_2_complex(1:min_f,1:min_t) - st_x_avg_2_complex(1:min_f,1:min_t); 
    S_abs = st_x_ind_2_abs(1:min_f,1:min_t) - st_x_avg_2_abs(1:min_f,1:min_t);
    if (out_type >= 3) S_log = st_x_ind_2_log(1:min_f,1:min_t) - st_x_avg_2_log(1:min_f,1:min_t); end;
else % ind
    S_complex = st_x_ind_2_complex(1:min_f,1:min_t);
    S_abs = st_x_ind_2_abs(1:min_f,1:min_t);
    if (out_type >= 3) S_log = st_x_ind_2_log(1:min_f,1:min_t); end;
end    

[inv_S, t_inv] = stockwell_inverse_transform(abs(S_abs), angle(S_complex), rate);


if (st_type == 1)
	
	try 
    	is_gmn = h1_struct.st_struct.is_gmn;

		if (h1_struct.st_struct.is_gmn == 1)
		
			channel_idx = get_channel_index(h1_struct.st_struct.channel_list, chan_name);
			S_ind_abs = h1_struct.data_struct.hdf1_ind_st_abs_S{channel_idx};
			
			if (apply_baseline == 1) 	
					
				if (out_type >= 3)
					log_transform = 1;
					S_ind_log = apply_tfr_abs_baseline(S_ind_abs, st_baseline_time_min_ms, st_baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type, log_transform);
				else
					log_transform = 0;
					S_ind_abs = apply_tfr_abs_baseline(S_ind_abs, st_baseline_time_min_ms, st_baseline_time_max_ms, rate, pre_stim_time_ms, baseline_type, log_transform);
				end	
			else
				if (out_type >= 3)
					S_ind_abs = log(abs(S_ind_abs));
				end			
			end
			
			S_abs = S_ind_abs(1:min_f,1:min_t);
			if (out_type >= 3) S_log = S_ind_log(1:min_f,1:min_t); end;
			
		end

	catch
		disp 'Warning: using mean inverse_ind_S';
	end
end


% calculate power from S amplitude

if (out_type == 1)
	if (apply_baseline == 1)
		S = (S_abs.^2) .* sign(S_abs);
	else
    	S = S_abs.^2;
	end
elseif (out_type == 2)
    S = S_abs;    
elseif (out_type == 3) % log(x^2) = 2*log(x)
    S = 2*S_log;	
elseif (out_type == 4)
    S = S_log;
else
    if (apply_baseline == 1)
		S = (S_abs.^2) .* sign(S_abs);
	else
    	S = S_abs.^2;
	end
end    

