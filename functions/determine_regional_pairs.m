function [pairs, pair_inds, pair_indlbls]=determine_regional_pairs(coords_filepath,relation,regions,sparse_logic)

% relation is a string, can be 'inter', 'intra', or 'both'
% regions is a cell, can contain strings 'frontal', 'central', 'parietal', ...
% 'occipital', 'right temporal', and/or 'left temporal'

%coords_filepath='/export/home/mike/matlab/origin/coords/61chans_ns.mat';
%relation='both';
%regions={'frontal','parietal','occipital'};
%regions={'frontal','central','parietal','occipital','lefttemporal','righttemporal'};


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

if n_channels ~= 61
    error('Regional coherence pairs available only for the 61 channel set.')
    return
end

if sparse_logic
    frontal=[5 6 7 8 9 34 35 43 44 47];
    central=[16 17 18 39 40 41 42 52 53 56];
    parietal=[21 22 23 24 25 48 49 59 60 61];
    occipital=[28 29 30 31 54 55 57 58];
    lefttemporal=[15 19 27 36 46 50];
    righttemporal=[14 20 26 37 45 51];
else
    frontal=[1 2 3 4 5 6 7 8 9 32 33 34 35 38 43 44 47];
    central=[10 11 12 13 16 17 18 39 40 41 42 52 53 56];
    parietal=[21 22 23 24 25 48 49 59 60 61];
    occipital=[28 29 30 31 54 55 57 58];
    lefttemporal=[15 19 27 36 46 50];
    righttemporal=[14 20 26 37 45 51];
end

region_names={'frontal','central','parietal','occipital','lefttemporal','righttemporal'};
region_chankey={frontal, central, parietal, occipital, lefttemporal, righttemporal};

pairs=[];
pair_inds=[];
pair_indlbls={};
relation_dummy=0;

n_regions=length(regions);
region_inds=zeros(n_regions,1);

for region=1:n_regions
    region_inds(region)=find(strcmpi(region_names,regions{region}));
end

%intra computes all pairs inside a region

if strcmpi(relation,'intra') || strcmpi(relation,'both')
for region=1:n_regions

    intra_pairs=nchoosek(region_chankey{region_inds(region)},2);
    
    relation_dummy=relation_dummy+1;
    intra_inds=ones(length(intra_pairs),1)*relation_dummy;

    pairs=[pairs;intra_pairs];
    pair_inds=[pair_inds; intra_inds];
    pair_indlbls{relation_dummy}=region_names{region_inds(region)};

end
end

%inter computes all pairs between two sets of nodes (complete bipartite
%graph)

if strcmpi(relation,'inter') || strcmpi(relation,'both')

biparts=nchoosek(region_inds,2);
n_biparts=size(biparts,1);

for bipart=1:n_biparts
    
    from_chans=region_chankey{biparts(bipart,1)};
    to_chans=region_chankey{biparts(bipart,2)};
    
    b_pairs=bipartition(from_chans,to_chans);
    
    relation_dummy=relation_dummy+1;
    b_inds=ones(length(b_pairs),1)*relation_dummy;
    
    pairs=[pairs;b_pairs];
    pair_inds=[pair_inds; b_inds];
    pair_indlbls{relation_dummy}=[  region_names{biparts(bipart,1)}, '-', ...
        region_names{biparts(bipart,2)} ] ;

end
end

end