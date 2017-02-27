function [h1_struct] = create_epoch_average_data(h1_struct)
%
% function [h1_struct] = create_epoch_average_data(h1_struct)
%

if (h1_struct.file_struct.data_type ~= 0)  
    return
end

h1_struct.data_struct.hdf1_avg_data = 0;

case_num = h1_struct.case_struct.case_num;
  
for the_case = 1:length(case_num)
    case_number = case_num(the_case);
    
    case_accepted_idx = get_accepted_case_trial_idxs(h1_struct, case_number);
    
    the_case_data = squeeze(mean(h1_struct.data_struct.hdf1_epoch_data(case_accepted_idx,:,:),1));  
    
    for channel=1:size(the_case_data,1)
       for sample=1:size(the_case_data,2)
           h1_struct.data_struct.hdf1_avg_data(the_case,channel,sample) = the_case_data(channel,sample);
       end    
    end
    % h1_struct.data_struct.hdf1_avg_data(the_case,:,:) = the_case_data(:,:);
end    

h1_struct.experiment_struct.n_samples = size(the_case_data,2);
