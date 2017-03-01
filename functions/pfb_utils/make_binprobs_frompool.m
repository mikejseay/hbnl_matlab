function bin_probs = make_binprobs_frompool(pool, n_bins)
% make an n_bins-long set of event probabilities from a pool of fixed
% values, pool

bin_probs = zeros(n_bins, 1);

n_elems = length(pool);
bins_perelem = n_bins ./ n_elems;
pool_inds = 1:n_bins;

for e = 1:n_elems

    pool_len = length(pool_inds);
    
    elem_inds = sort(randperm(pool_len, bins_perelem));
    elem_inds = pool_inds(elem_inds);
    
    bin_probs(elem_inds) = pool(e);
    pool_inds = setdiff(pool_inds, elem_inds);
    
end

end