function [S_values, abs_S, inv_S, t, f] = calc_chan_stockwell_values(h1_struct, min_frequencies, max_frequencies, min_time_ms, max_time_ms, value_type, chan_name, st_type, apply_baseline, baseline_type, st_baseline_time_min_ms, st_baseline_time_max_ms, out_type)
%
% function [S_values, S, inv_S, t, f] = calc_chan_stockwell_values(h1_struct,
% min_frequencies, max_frequencies, min_time_ms, max_time_ms, value_type, 
% chan_name, [st_type], [apply_baseline], [baseline_type], [st_baseline_time_min_ms], 
% [st_baseline_time_max_ms], [out_type])
%
if (nargin < 13)
   out_type = 1;
end 
if (nargin < 12)
   out_type = 1;
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
end 
if (nargin < 11)
   out_type = 1; 
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
   st_baseline_time_min_ms = h1_struct.st_struct.st_baseline_time_min_ms;
end 
if (nargin < 10)
   out_type = 1; 
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
   st_baseline_time_min_ms = h1_struct.st_struct.st_baseline_time_min_ms;
   baseline_type = h1_struct.st_struct.baseline_type;
end 
if (nargin < 9)
   out_type = 1; 
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
   st_baseline_time_min_ms = h1_struct.st_struct.st_baseline_time_min_ms;
   baseline_type = h1_struct.st_struct.baseline_type;
   apply_baseline = 0;
end 
if (nargin < 8)
   out_type = 1; 
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
   st_baseline_time_min_ms = h1_struct.st_struct.st_baseline_time_min_ms;
   baseline_type = h1_struct.st_struct.baseline_type;
   apply_baseline = 0;
   st_type = 0;
end  

rate = h1_struct.experiment_struct.rate;

pre_stim_time_secs = h1_struct.experiment_struct.pre_stim_time_ms / 1000;

[abs_S, inv_S, t, f] = calc_chan_stockwell_power(h1_struct, chan_name, st_type, apply_baseline, baseline_type, st_baseline_time_min_ms, st_baseline_time_max_ms, out_type);

n_freqs = min(min([length(min_frequencies) length(min_frequencies)]));

min_time_secs = (min_time_ms / 1000) + pre_stim_time_secs;
max_time_secs = (max_time_ms / 1000) + pre_stim_time_secs;

if (min_time_secs > max_time_secs)
   min_time_secs = t(1);
   max_time_secs = t(length(t));
end   
if (min_time_secs <= 0)
   min_time_secs = t(1);
end
if (max_time_secs <= 0)
   max_time_secs = t(length(t));
end
if (max_time_secs > t(length(t)))
   max_time_secs = t(length(t));
end

[junk, min_time_idx] = min(abs(t - min_time_secs));
[junk, max_time_idx] = min(abs(t - max_time_secs));

for i=1:n_freqs

    min_frequency = min_frequencies(i);
    max_frequency = max_frequencies(i);
    
    if (min_frequency > max_frequency)
       min_frequency = f(2);
       max_frequency = f(length(f));
    end   
    if (min_frequency <= 0)
       min_frequency = f(2);
    end
    if (max_frequency <= 0)
       max_frequency = f(length(f));
    end
    if (max_frequency > f(length(f)))
       max_frequency = f(length(f));
    end

    [junk, min_frequency_idx]= min(abs(f - min_frequency));
    [junk, max_frequency_idx] = min(abs(f - max_frequency));

    if (min_frequency_idx == 1)
       min_frequency_idx = 2; 
    end 
    if (f(min_frequency_idx) == 0.0)
       min_frequency_idx = min_frequency_idx + 1; 
    end  
    
    if (max_frequency_idx > length(f))
       max_frequency_idx = f(length(f)); 
    end 
   
    edge_length_ms = 100.0;
    abs_S = mask_tfr_matrix_edges(abs_S, edge_length_ms, rate);
    
    if (value_type == 1)
        S_values{i} = mean(mean(abs_S(min_frequency_idx:max_frequency_idx,min_time_idx:max_time_idx)));    
    elseif (value_type == 2)
        S_values{i} = max(max(abs_S(min_frequency_idx:max_frequency_idx,min_time_idx:max_time_idx)));
    elseif (value_type == 3)
        S_values{i} = get_matrix_centroid_value(abs_S(min_frequency_idx:max_frequency_idx,min_time_idx:max_time_idx));
    elseif (value_type == 4)
        [junk, freq_idx] = max(max(abs_S(min_frequency_idx:max_frequency_idx,min_time_idx:max_time_idx),[],2));
        freq_idx = freq_idx + min_frequency_idx; 
		S_values{i} = f(freq_idx);
    elseif (value_type == 5)
        [junk, time_idx] = max(max(abs_S(min_frequency_idx:max_frequency_idx,min_time_idx:max_time_idx),[],1));
		time_idx = time_idx + min_time_idx;
		S_values{i} = (t(time_idx) - pre_stim_time_secs) * 1000.0;  
    elseif (value_type == 6)
        S_values{i} = mean(mean(abs_S(min_frequency_idx:max_frequency_idx,min_time_idx:max_time_idx)));    
    else    
        S_values{i} = mean(mean(abs_S(min_frequency_idx:max_frequency_idx,min_time_idx:max_time_idx)));
    end
    
end

