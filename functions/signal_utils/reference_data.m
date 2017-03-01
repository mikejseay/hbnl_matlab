function [ data_out ] = reference_data( data_in , scheme )
%REFERENCE_DATA References data according to an input parameter
% Inputs should be data that is n_samples x n_channels, and a scheme should
% be specified. 'average' or 'monopolar' are accepted.

    [n_samps,n_chans]=size(data_in);
    
    %average referencing
    if strcmpi(scheme,'average')
        %first check if there are any extreme channels
        skip_pts=500; %number of edge points to exclude from examination
        mean_thresh=10; % threshold for channel mean amplitude
        std_thresh=500; % threshold for channel standard deviation
        bad_by_mean=find(mean(data_in(1+skip_pts:end-skip_pts,:),1)>mean_thresh);
        bad_by_std=find(std(data_in(1+skip_pts:end-skip_pts,:),1)>std_thresh);
        bad_channels=union(bad_by_mean,bad_by_std);
        good_channels=setdiff(1:n_chans,bad_channels);
        %then reference to the average of all non-extreme electrodes
       data_out = bsxfun(@minus, data_in, mean(data_in(:,good_channels),2));
   %no referencing
    else
        data_out=data_in;
    end
