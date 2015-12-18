function plot_pairs(pairmat,chan_locs,pair_inds,clickable,direction_spec)
% plot a set of pairs (n_pairs x 2) as arcs on a 2d headplot

% inputs
% ------
% pairmat: 2-column array (n_pairs x 2), #'s refer to channel #'s
% chan_locs: channel locations structure in eeglab style, #'s match pairs above
% pair_inds: numbers that group pairs into subsets
% clickable: logical, whether the plot will have an interative clickable legend
% direction_spec: 1, 2, or 3. aesthetics of arc plots

% written by michael seay, hbnl, 2015

if nargin<3
    pair_inds=[];
end

if nargin<4
    clickable=false;
end

if nargin<5
    direction_spec=1;
end

n_chosen_pairs=1:size(pairmat,1);

if isempty(pair_inds)
    hyp_colors=false;
    distant_pair_colors=distinguishable_colors(length(n_chosen_pairs));
else
    % make colors based on pair hypothesis
    hyp_colors=true;
    n_pairtypes=max(pair_inds);
    distant_pair_colors=distinguishable_colors(n_pairtypes);
end

%plot them as arcs on a head
dummy_data=zeros(length(chan_locs),1);
topoplot(dummy_data,chan_locs,'style','blank'); hold on;
for chosen_pair=n_chosen_pairs
    %define [x1 x2], [y1 y2]
    x=[chan_locs(pairmat(chosen_pair,1)).topo_x chan_locs(pairmat(chosen_pair,2)).topo_x];
    y=[chan_locs(pairmat(chosen_pair,1)).topo_y chan_locs(pairmat(chosen_pair,2)).topo_y];
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
    if hyp_colors
        arc_h(chosen_pair)=Draw_Arc_Clockwise([x(1),y(1)], [x(2),y(2)], ...
            distant_pair_colors(pair_inds(chosen_pair),:), 1, direction); hold on;
    else
        arc_h(chosen_pair)=Draw_Arc_Clockwise([x(1),y(1)], [x(2),y(2)], ...
            distant_pair_colors(chosen_pair,:), 1, direction); hold on;
    end
end

%make labels
p_label=cell(length(n_chosen_pairs),1);
for pair=n_chosen_pairs
    p_label{pair}=[chan_locs(pairmat(pair,1)).labels,'-',chan_locs(pairmat(pair,2)).labels];
end

if clickable
    clickableLegend(arc_h,{p_label{n_chosen_pairs}})
end

end