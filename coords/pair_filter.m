function [outpairs,outinds]=pair_filter(inpairs,ininds,chan_locs,filter_type,arg)
% filters an array of pairs in one of several modes

% inputs
% ------
% inpairs: array of pairs (n_pairs x 2)
% ininds: array of grouping indices corresponding to inpairs (n_pairs x 1)
% chan_locs: channel locations structure in eeglab style, #'s match pairs above
% filter_type: type of filter being used
% arg: argument for filter, see guide below

% outputs
% -------
% outpairs: filtered pairs
% outinds: filtered grouping indices

% filter types
% ------------
% 'x_edge' excludes outer-most electrodes (on edges of head in 2d headplot)
    % no arg required
    
% 'distance' filters based on distance between channels in each pair
    % arg should be [lower upper] in cm OR
    % a single number,
        % 1 = < mean - std
        % 2 = (mean - std , mean)
        % 3 = (mean - std, mean + std)
        % 4 = (mean, mean + std)
        % 5 = > mean + std
        
% 'angle' filters based on orientation of electrode pairs
    % arg should be:
    % 'ap' = anterior-posterior oriented pairs
    % 'lr' = left-right oriented pairs
    % 'aplr' = vertical and horizontal pairs (nearly anterior-posterior or
    %left-right)
    % 'diag' = diagonal pairs
    % 'all' = all pairs regardless of orientation (don't filter)
    
% 'max_degree' filters pairs such that each electrode is connected to at most
% arg other electrodes
    % arg should be 1 number
    
% 'hemispheric' filters pairs based on hemispheric relationship
    % arg should be
    % 'intra' = only intrahemispheric pairs (e.g. F1-F3)
    % 'inter' = only interhemispheric pairs (e.g. F3-F4)
    % 'neither' = only pairs that contain a midline electrode (e.g. FZ-F5)
    % 'all' = return all (don't filter)

% written by michael seay, hbnl, 2015

n_pairs=size(inpairs,1);
n_chans=length(chan_locs);

switch filter_type
    case 'x_edge'
        
        exclude_channels=[1 2 3 4 14 15 26 27 30 31 32 33 36 37 38 46 45 54 55 58];
        exclude_pairs=any(ismember(inpairs,exclude_channels),2);
        outpairs=inpairs(~exclude_pairs,:);
        outinds=ininds(~exclude_pairs);
        
    case 'distance'
        
        if ~isnumeric(arg)
            error('For distance filter, arg must be [lower upper] in cm')
            return
        end
        
        p_dists=zeros(n_pairs,1);
        for p=1:n_pairs
            x_dist = chan_locs(inpairs(p,1)).X - chan_locs(inpairs(p,2)).X;
            y_dist = chan_locs(inpairs(p,1)).Y - chan_locs(inpairs(p,2)).Y;
            z_dist = chan_locs(inpairs(p,1)).Z - chan_locs(inpairs(p,2)).Z;
            p_dists(p) = sqrt( x_dist^2 + y_dist^2 + z_dist^2 );
        end
        
        if length(arg)==2
        
        arg = arg / 10; %convert from cm to dm
        distance_inds = p_dists > arg(1) & p_dists < arg(2);
        
        elseif length(arg)==1
        
        %determine subsets by distance
        mean_distance = mean(p_dists);
        std_distance = std(p_dists);

        switch arg
            case 1 % 1 = very near
                distance_inds = p_dists < mean_distance - std_distance;
            case 2 % 2 = lower middle distance (mean-std, mean)
                distance_inds = p_dists > mean_distance - std_distance & p_dists < mean_distance;
            case 3 % 2 = middle distance (mean-std, mean+std)
                distance_inds = p_dists < mean_distance + std_distance & p_dists > mean_distance - std_distance;
            case 4 % 3 = upper middle distance (mean,mean+std)
                distance_inds = p_dists > mean_distance & p_dists < mean_distance + std_distance;
            case 5 % 4 = very distant
                distance_inds = p_dists > mean_distance + std_distance;
        end
            
        end
        
        outpairs=inpairs(distance_inds,:);
        outinds=ininds(distance_inds);
        
    case 'angle'
        
        %determine the angle of the chosen pairs
        pair_angle=zeros(n_pairs,1);
        for p=1:n_pairs
            %define [x1 x2], [y1 y2]
            x=[chan_locs(inpairs(p,1)).topo_x chan_locs(inpairs(p,2)).topo_x];
            y=[chan_locs(inpairs(p,1)).topo_y chan_locs(inpairs(p,2)).topo_y];
            pair_angle(p)=atan((y(2)-y(1))/(x(2)-x(1)));
        end
        
        oriented_inds = pair_angle < -1.26 | ( pair_angle > -0.314 & pair_angle < 0.314 ) | pair_angle > 1.26 ;
        ap_inds = pair_angle < -1.26 | pair_angle > 1.26 ;
        lr_inds = pair_angle > -0.314 & pair_angle < 0.314;
        %take the subset of non-nearly a-p or l-r oriented pairs
        non_oriented_inds = ~oriented_inds;
        switch arg
            case 'ap'
                angle_inds = ap_inds;
            case 'lr'
                angle_inds = lr_inds;
            case 'aplr' %take the a-p or l-r oriented pairs
                angle_inds = oriented_inds;
            case 'diag' %take the diagonally-oriented pairs
                angle_inds = non_oriented_inds;
            case 'all'
                angle_inds = 1:length(pair_angle);
        end
        outpairs = inpairs(angle_inds,:);
        outinds = ininds(angle_inds);
        
    case 'max_degree'
        
        unique_chans=unique(inpairs)';
        degree=zeros(n_chans,1);
        
        % this algorithm removes pairs until the degree is reduced to a certain
        % number. we will sort the indices by distance, and remove
        % every nth pair of that channel membership to get "even" amounts of
        % pair removal with distance
        
        % honestly this algorithm is pretty janky but it works to reduce
        % the total number of pairs
        
        p_dists=pair_distance(inpairs,chan_locs);
        
        [~,sort_inds]=sort(p_dists);
        unsort_inds(sort_inds)=1:length(p_dists);
        
        sorted_pairs=inpairs(sort_inds,:);
        remlogic=false(size(inpairs,1),n_chans);
        
        for chan=unique_chans
            [memberpair_inds,~]=find(any(sorted_pairs==chan,2));
            n_memberpairs=length(memberpair_inds);
            n_pairs2rem = n_memberpairs - arg;
            step=floor(n_memberpairs/n_pairs2rem);
            if ~isinf(step)
                rem_inds=memberpair_inds(1:step:end);
                remlogic(rem_inds(1:n_pairs2rem),chan)=true;
            end
        end
        
        pair_rem=any(remlogic,2);
        pair_rem=pair_rem(unsort_inds);
        outpairs=inpairs(~pair_rem,:);
        outinds=ininds(~pair_rem);
        
        
    case 'hemispheric'
        
        hemi_relation=zeros(n_pairs,1);
        % if y coord is positive, left, if y coord is negative, right
        for p=1:n_pairs
            if ( chan_locs(inpairs(p,1)).Y > 0 && chan_locs(inpairs(p,2)).Y > 0 ) || ...
                    ( chan_locs(inpairs(p,1)).Y < 0 && chan_locs(inpairs(p,2)).Y < 0 ) % if pairs have the same sign (both on left, both on right)
                hemi_relation(p) = 1; % 1 means intra hemispheric
            elseif chan_locs(inpairs(p,1)).Y==0 || chan_locs(inpairs(p,2)).Y==0 % if either pair is at 0 (on midline)
                hemi_relation(p) = 3; % 3 means neither, or relating to the midline
            else
                hemi_relation(p) = 2; % 2 means inter hemispheric
            end
        end
        
        switch arg
            case 'intra'
                pair_inds = hemi_relation == 1;
            case 'inter'
                pair_inds = hemi_relation == 2;
            case 'neither'
                pair_inds = hemi_relation == 3;
            case 'all'
                pair_inds = 1:n_pairs;
        end
        
        outpairs=inpairs(pair_inds,:);
        outinds=ininds(pair_inds);
        
end
    
    
end