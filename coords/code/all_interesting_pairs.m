pairs = nchoosek(1:61,2);
pair_inds = ones(size(pairs));
load('/export/home/mike/matlab/origin/coords/61chans_ns.mat')
[pairs2,pair_inds2] = pair_filter(pairs,pair_inds,chan_locs,'distance',[5 Inf]);
[pairs3,pair_inds3] = pair_filter(pairs2,pair_inds2,chan_locs,'angle','aplr');
[pairs4,pair_inds4] = pair_filter(pairs3,pair_inds3,chan_locs,'max_degree',12);
plot_pairs(pairs4,chan_locs);