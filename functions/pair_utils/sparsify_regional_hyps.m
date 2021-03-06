function [outpairs, outinds] = sparsify_regional_hyps( inpairs, ininds, ...
    chan_locs, maxpairs_perhyp, start_degree )

if nargin < 5
    start_degree = 18;
end

n_hyps = max(ininds);

outpairs = [];
outinds = [];
for hyp = 1:n_hyps
p = inpairs(ininds==hyp, :);
n_p = size(p, 1);
p_i = hyp * ones(n_p, 1);
if n_p > maxpairs_perhyp
    for degree = start_degree:-1:1
        [p2, p_i2] = pair_filter(p, p_i, chan_locs, 'max_degree', degree);
        if size(p2, 1) <= maxpairs_perhyp; break; end
    end
    outpairs = [outpairs; p2];
    outinds = [outinds; p_i2];
else    
    outpairs = [outpairs; p];
    outinds = [outinds; p_i];
end

end