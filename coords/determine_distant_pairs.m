function pair_subset=determine_distant_pairs(coords_filepath,exclude_logic, ...
    distance_spec,orientation_spec,edges_per_chan,do_plot)

% determine pairs of electrodes for coherence analyses.

    %coords_filepath is number of electrodes (31 or 61) or the filepath of
    %a locations file. this has to be a chan_locs structure in the EEGLAB
    %style
    
    %exclude_logic indicates to exclude edge electrodes (1) or not (0)
    
    %distance_spec indicates distance of electrodes to examine in increasing
    % manner like so:
        % 1 = < mean - std
        % 2 = (mean - std , mean)
        % 3 = (mean - std, mean + std)
        % 4 = (mean, mean + std)
        % 5 = > mean + std
        
    %orientation_spec indicates orientation of electrodes to examine like
    %so:
        % 'all' = all pairs regardless of orientation
        % 'ap' = anterior-posterior oriented pairs
        % 'lr' = left-right oriented pairs
        % 'aplr' = vertical and horizontal pairs (nearly anterior-posterior or
        %left-right)
        % 'diag' = diagonal pairs
        
    %edges_per_chan indicates how many edges you desire per channel
        
    %direction_spec indicates how to plot the arcs (purely aesthetic)
        % 1 = all clockwise
        % 2 = alternate every other
        % 3 = alternate erratically
        
    % output is the indices of the chosen / plotted pairs
    % as a pairs x 2 matrix, can be used as the input for erp_analysis in
    % the param_struct.

%%
    
if isscalar(coords_filepath)
    %if coords spec is a number, load the file with that number of channels
    switch coords_filepath
        case 31
            load('/export/home/mike/matlab/origin/coords/31chans_ns.mat')
        case 61
            load('/export/home/mike/matlab/origin/coords/61chans_ns.mat')
    end
elseif ischar(coords_filepath)
    %if a string, load that location file
    load(coords_filepath)
elseif isstruct(coords_filepath)
    %if a structure, set that as the chan_locs
    chan_locs=coords_filepath;
else
    error('Channel locations incorrectly specified.')
end

n_channels=length(chan_locs);

%specify sub-set of electrodes to exclude as pairs (edge electrodes)
if exclude_logic
    if n_channels==31
        exclude_channels=[1 2 3 4 14 15 26 27 30 31];
    elseif n_channels==61
        exclude_channels=[1 2 3 4 14 15 26 27 30 31 32 33 36 37 38 46 45 54 55 58];
    else
        error('Channel locations file specified has unexpected # of electrodes.')
    end
else
    exclude_channels=[];
end

channel_subset=setdiff(1:n_channels,exclude_channels);

%determine pair distances
pairs=nchoosek(channel_subset,2);
n_pairs=size(pairs,1);
pair_distance=zeros(n_pairs,1);
for pair=1:n_pairs
    x_dist = chan_locs(pairs(pair,1)).X - chan_locs(pairs(pair,2)).X;
    y_dist = chan_locs(pairs(pair,1)).Y - chan_locs(pairs(pair,2)).Y;
    z_dist = chan_locs(pairs(pair,1)).Z - chan_locs(pairs(pair,2)).Z;
    pair_distance(pair) = sqrt( x_dist^2 + y_dist^2 + z_dist^2 );
end

if length(distance_spec)==1

%determine subsets by distance
mean_distance = mean(pair_distance);
std_distance = std(pair_distance);

switch distance_spec
    case 1 % 1 = very near
        distance_inds = pair_distance < mean_distance - std_distance;
    case 2 % 2 = lower middle distance (mean-std, mean)
        distance_inds = pair_distance > mean_distance - std_distance & pair_distance < mean_distance;
    case 3 % 2 = middle distance (mean-std, mean+std)
        distance_inds = pair_distance < mean_distance + std_distance & pair_distance > mean_distance - std_distance;
    case 4 % 3 = upper middle distance (mean,mean+std)
        distance_inds = pair_distance > mean_distance & pair_distance < mean_distance + std_distance;
    case 5 % 4 = very distant
        distance_inds = pair_distance > mean_distance + std_distance;
end

elseif length(distance_spec)==2
    distance_spec = distance_spec / 10; %convert from cm to dm
    distance_inds = pair_distance > distance_spec(1) & pair_distance < distance_spec(2);
end

%take distance subset
pair_subset = pairs(distance_inds,:);
n_chosen_pairs = size(pair_subset,1);
%distance_range = [min(pair_distance(distance_inds)) max(pair_distance(distance_inds))];

%determine the angle of the chosen pairs
distant_pair_angle=zeros(n_chosen_pairs,1);
for chosen_pair=1:n_chosen_pairs
    %define [x1 x2], [y1 y2]
    x=[chan_locs(pair_subset(chosen_pair,1)).topo_x chan_locs(pair_subset(chosen_pair,2)).topo_x];
    y=[chan_locs(pair_subset(chosen_pair,1)).topo_y chan_locs(pair_subset(chosen_pair,2)).topo_y];
    distant_pair_angle(chosen_pair)=atan((y(2)-y(1))/(x(2)-x(1)));
end

%take the subset of nearly a-p oriented pairs
oriented_inds = distant_pair_angle < -1.26 | ( distant_pair_angle > -0.314 & distant_pair_angle < 0.314 ) | distant_pair_angle > 1.26 ;
ap_inds = distant_pair_angle < -1.26 | distant_pair_angle > 1.26 ;
lr_inds = distant_pair_angle > -0.314 & distant_pair_angle < 0.314;
%take the subset of non-nearly a-p or l-r oriented pairs
non_oriented_inds = ~oriented_inds;
switch orientation_spec
    case 'all' %do nothing, take all pairs regardless of orientation
    case 'ap'
        pair_subset = pair_subset(ap_inds,:);
    case 'lr'
        pair_subset = pair_subset(lr_inds,:);
    case 'aplr' %take the a-p or l-r oriented pairs
        pair_subset = pair_subset(oriented_inds,:);
    case 'diag' %take the diagonally-oriented pairs
        pair_subset = pair_subset(non_oriented_inds,:);
end

n_chosen_pairs=size(pair_subset,1);

%%

if edges_per_chan > 0
%sort pairs by distance
%determine the degree (number of connections per channel) and reduce to have equal
%connections across channels
unique_chans=unique(pair_subset)';
degree=zeros(n_channels,1);
%this algorithm removes pairs until the degree is reduced to a certain
%number
for chan=unique_chans
pair=0;
while pair<=length(pair_subset)-1
    pair=pair+1;
    if any(pair_subset(pair,:)==chan)
        degree(chan)=degree(chan)+1;
    end
    if degree(chan)>edges_per_chan
        if any(pair_subset(pair,:)==chan)
            pair_subset(pair,:)=[];
            pair=pair-1;
        end
    end
end
end

end

%%

if do_plot

plot_pairs(pair_subset,chan_locs,1);

end

fprintf('%d pairs determined.\n',size(pair_subset,1))

end