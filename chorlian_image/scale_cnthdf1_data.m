function [h1_struct] = scale_cnthdf1_data(h1_struct)
%
% [hdf1_struct] = scale_cnthdf1_data(hdf1_struct)
%
if (h1_struct.file_struct.data_type ~= 0)  
    return
end

h1_struct.data_struct.hdf1_cnt_data = ...
	h1_struct.data_struct.hdf1_cnt_data * h1_struct.run_struct.scale_value(1);
% for i=1:size(h1_struct.data_struct.hdf1_cnt_data,1)
%     h1_struct.data_struct.hdf1_cnt_data(i,:) = h1_struct.data_struct.hdf1_cnt_data(i,:) .* h1_struct.run_struct.scale_value(i);
% end
