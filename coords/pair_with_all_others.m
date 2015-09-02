%load location file
coords_filepath='A:\matlab\coords\31chans_ns.mat';
load(coords_filepath)
n_channels=length(chan_locs);

%specify electrodes to pair with all others
% FZ, CZ, PZ
electrode_seeds=[7,16,25];
n_seeds=length(electrode_seeds);

%generate pairs
pairs=zeros(n_seeds*(n_channels-1),2);
for seed=1:length(electrode_seeds)
    pairs( ((n_channels-1)*(seed-1)+1):((n_channels-1)*seed), : ) = ...
        [ ones(n_channels-1,1)*electrode_seeds(seed) , setdiff(1:n_channels,electrode_seeds(seed))' ];
end

%determine pair distances
n_pairs=size(pairs,1);
pair_distance=zeros(n_pairs,1);
for pair=1:n_pairs
    x_dist = chan_locs(pairs(pair,1)).X - chan_locs(pairs(pair,2)).X;
    y_dist = chan_locs(pairs(pair,1)).Y - chan_locs(pairs(pair,2)).Y;
    z_dist = chan_locs(pairs(pair,1)).Z - chan_locs(pairs(pair,2)).Z;
    pair_distance(pair) = sqrt( x_dist^2 + y_dist^2 + z_dist^2 );
end

%%

%determine subsets by distance
mean_distance = mean(pair_distance);
std_distance = std(pair_distance);
distant_inds = pair_distance > mean_distance + std_distance;
middle_inds = pair_distance < mean_distance + std_distance & pair_distance > mean_distance - std_distance;
upper_middle_inds = pair_distance < mean_distance + std_distance & pair_distance > mean_distance;
near_inds = pair_distance < mean_distance - std_distance;

%take a subset (don't because seeded)
pair_subset = pairs(:,:);
n_chosen_pairs = size(pair_subset,1);
distance_range = [min(pair_distance(upper_middle_inds)) max(pair_distance(upper_middle_inds))];

%determine the angle of the chosen pairs
distant_pair_angle=zeros(n_chosen_pairs,1);
for chosen_pair=1:n_chosen_pairs
    %define [x1 x2], [y1 y2]
    x=[chan_locs(pair_subset(chosen_pair,1)).topo_x chan_locs(pair_subset(chosen_pair,2)).topo_x];
    y=[chan_locs(pair_subset(chosen_pair,1)).topo_y chan_locs(pair_subset(chosen_pair,2)).topo_y];
    distant_pair_angle(chosen_pair)=atan((y(2)-y(1))/(x(2)-x(1)));
end

if false

    %take the subset of nearly a-p or l-r oriented pairs
    oriented_inds = distant_pair_angle < -1.26 | ( distant_pair_angle > -0.314 & distant_pair_angle < 0.314 ) | distant_pair_angle > 1.26 ;
    %take the subset of non-nearly a-p or l-r oriented pairs
    non_oriented_inds = ~oriented_inds;
    pair_subset = pair_subset(oriented_inds,:);
    n_chosen_pairs=size(pair_subset,1);

end

%assign them colors??
distant_pair_colors=distinguishable_colors(n_chosen_pairs);

%%

%plot them as arcs on a head
dummy_data=zeros(length(chan_locs),1);
for seed=1:n_seeds
    figure; topoplot(dummy_data,chan_locs,'style','blank'); hold on;
    for chosen_pair=((n_channels-1)*(seed-1)+1):((n_channels-1)*seed)
        %define [x1 x2], [y1 y2]
        x=[chan_locs(pair_subset(chosen_pair,1)).topo_x chan_locs(pair_subset(chosen_pair,2)).topo_x];
        y=[chan_locs(pair_subset(chosen_pair,1)).topo_y chan_locs(pair_subset(chosen_pair,2)).topo_y];
        %determine the direction
        %direction=mod(mod(chosen_pair,3),2);
        %direction=mod(chosen_pair,2);
        direction=1;
        %plot the arc
        Draw_Arc_Clockwise([x(1),y(1)], [x(2),y(2)], distant_pair_colors(chosen_pair,:), 1, direction); hold on;
        %Draw_Arc_Clockwise([x(1),y(1)], [x(2),y(2)], [0 0 0], 1, direction); hold on;
    end
end