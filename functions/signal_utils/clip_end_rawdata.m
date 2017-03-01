function [clipped_data] = clip_end_rawdata(data,thresh)

    n_samps = size(data,1);
    for samp = n_samps:-1:1
        if std(data(samp,:))<thresh
            clipped_data=data(1:samp,:);
            return
        end
    end
    clipped_data = data;
    
