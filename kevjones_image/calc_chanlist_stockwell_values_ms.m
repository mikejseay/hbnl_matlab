function [S_chan_values, Sbase_chan_values, S_all, t, f] = calc_chanlist_stockwell_values_ms(h1_struct, ...
    min_frequencies, max_frequencies, min_time_ms, max_time_ms, value_type, chan_name_list, st_type, ...
    apply_baseline, baseline_type, st_baseline_time_min_ms, st_baseline_time_max_ms, out_type)
%
% function [S_chan_values, Sbase_chan_values] = calc_chanlist_stockwell_values(h1_struct, min_frequencies, 
% max_frequencies, min_time_ms, max_time_ms, value_type, [chan_name_list], [st_type], 
% [apply_baseline], [baseline_type], [st_baseline_time_min_ms],
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
if (nargin < 7)
   out_type = 1; 
   st_baseline_time_max_ms = h1_struct.st_struct.st_baseline_time_max_ms;
   st_baseline_time_min_ms = h1_struct.st_struct.st_baseline_time_min_ms;
   baseline_type = h1_struct.st_struct.baseline_type;
   apply_baseline = 0;
   st_type = 0;
   chan_name_list = h1_struct.st_struct.channel_list;
end  

n_chans = size(h1_struct.st_struct.channel_list,2);

%%% NEW CODE FOR SAVING 3D ARRAYS OF S_VALUES %%%
tmp_mean_data = zeros(1, size(h1_struct.data_struct.hdf1_avg_st_data, 2));
rate = h1_struct.experiment_struct.rate;
[tmp_S,tmp_t,tmp_f] = stockwell_transform(tmp_mean_data, 0, length(tmp_mean_data)/2, 1/rate);
clear tmp_mean_data tmp_f tmp_t rate
S_all = zeros(size(tmp_S, 1), size(tmp_S, 2), n_chans);
%%% NEW CODE FOR SAVING 3D ARRAYS OF S_VALUES %%%
for j=1:n_chans
    
    chan_name = upper(strtrim(h1_struct.st_struct.channel_list{j}));
    
    [S_values, S, inv_S, t, f] = calc_chan_stockwell_values_ms(h1_struct, min_frequencies, max_frequencies, min_time_ms, max_time_ms, value_type, chan_name, st_type, apply_baseline, baseline_type, st_baseline_time_min_ms, st_baseline_time_max_ms, out_type);  
    % [Sbase_values, S, inv_S, t, f] = calc_chan_stockwell_values(h1_struct, min_frequencies, max_frequencies, st_baseline_time_min_ms, st_baseline_time_max_ms, value_type, chan_name, st_type, 0, baseline_type, st_baseline_time_min_ms, st_baseline_time_max_ms, out_type);
    Sbase_values = 0;
    
    S_chan_values{j} = S_values;
    Sbase_chan_values{j} = Sbase_values;
    
    %%% NEW CODE FOR SAVING 3D ARRAYS OF S_VALUES %%%
    S_all(:, :, j) = S;
    %%% NEW CODE FOR SAVING 3D ARRAYS OF S_VALUES %%%
end


