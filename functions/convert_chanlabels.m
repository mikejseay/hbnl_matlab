function inds = convert_chanlabels(chanlocs, labels)
% convert channel label (strings) into indices

inds=[];
for c=1:length(labels)
    tmp_ind = find(strcmpi({chanlocs.labels}, labels{c}));
    inds=[inds,tmp_ind];
end

end