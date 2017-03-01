% create an eeglab-compatiable channel locations data structure with only a
% subset of the total 61 neuroscan channels

load('/export/home/mike/matlab/coords/61chans_ns.mat')
good_inds=[4 5 6 7 8 11 12 15 39 40 43 44 47 52 53 56];
good_inds(good_inds<31)=good_inds(good_inds<31)+1;
bad_inds=setdiff(1:61,good_inds);
bad_inds=fliplr(bad_inds);
chan_locs(bad_inds)=[];