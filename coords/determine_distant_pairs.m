function pair_subset=determine_distant_pairs(coords_filepath,exclude_logic, ...
    distance_spec,orientation_spec,edges_per_chan,direction_spec)

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
        % 1 = all pairs regardless of orientation
        % 2 = vertical/horizontal pairs (nearly anterior-posterior or
        %left-right)
        % 3 = diagonal pairs
        
    %direction_spec indicates how to plot the arcs (purely aesthetic)
        % 1 = all clockwise
        % 2 = alternate every other
        % 3 = alternate erratically
        
    %edges_per_chan indicates how many edges you desire per channel
        
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
%take the subset of non-nearly a-p or l-r oriented pairs
non_oriented_inds = ~oriented_inds;
switch orientation_spec
    case 1 %do nothing, take all pairs regardless of orientation
    case 2 %take the a-p or l-r oriented pairs
        pair_subset = pair_subset(oriented_inds,:);
    case 3 %take the diagonally-oriented pairs
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
n_chosen_pairs=1:size(pair_subset,1);
%assign them colors
%n_chosen_pairs=31:60;
distant_pair_colors=distinguishable_colors(length(n_chosen_pairs));
%n_chosen_pairs=31:60;

%plot them as arcs on a head
dummy_data=zeros(length(chan_locs),1);
figure; topoplot(dummy_data,chan_locs,'style','blank'); hold on;
for chosen_pair=n_chosen_pairs
    %define [x1 x2], [y1 y2]
    x=[chan_locs(pair_subset(chosen_pair,1)).topo_x chan_locs(pair_subset(chosen_pair,2)).topo_x];
    y=[chan_locs(pair_subset(chosen_pair,1)).topo_y chan_locs(pair_subset(chosen_pair,2)).topo_y];
    %determine the direction
    switch direction_spec
        case 1
            direction=1;
        case 2
            direction=mod(chosen_pair,2);
        case 3
            direction=mod(mod(chosen_pair,3),2);
    end
    %plot the arc
    arc_h(chosen_pair)=Draw_Arc_Clockwise([x(1),y(1)], [x(2),y(2)], distant_pair_colors(chosen_pair,:), 3, direction); hold on;
end

%make labels
p_label=cell(length(n_chosen_pairs),1);
for pair=n_chosen_pairs
    p_label{pair}=[chan_locs(pair_subset(pair,1)).labels,'-',chan_locs(pair_subset(pair,2)).labels];
end

clickableLegend(arc_h,{p_label{n_chosen_pairs}})

fprintf('%d pairs determined.\n',size(pair_subset,1))

end