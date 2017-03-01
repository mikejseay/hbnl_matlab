function distances=pair_distance(inpairs,chan_locs)
% calculate distance for each pair of channels in array (n_pairs x 2)

% inputs
% ------
% inpairs: n_pairs x 2 array of channel indices for each pair
% chan_locs: channel locations structure in eeglab style, #'s match pairs above

% written by michael seay, hbnl, 2015

n_pairs = size(inpairs,1);
distances = zeros(n_pairs,1);
for p=1:n_pairs
    x_dist = chan_locs(inpairs(p,1)).X - chan_locs(inpairs(p,2)).X;
    y_dist = chan_locs(inpairs(p,1)).Y - chan_locs(inpairs(p,2)).Y;
    z_dist = chan_locs(inpairs(p,1)).Z - chan_locs(inpairs(p,2)).Z;
    distances(p) = sqrt( x_dist^2 + y_dist^2 + z_dist^2 ) * 10; %convert to cm
end

end