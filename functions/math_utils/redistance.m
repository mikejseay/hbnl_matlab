function p_dists_jit = redistance(p_dists)

[p_dists_sorted,sort_inds]=sort(p_dists);
unsort_inds(sort_inds)=1:length(p_dists);
p_dists_jit = unsort_inds';
p_dists_jit = norm2limits_mat(p_dists_jit, [min(p_dists) max(p_dists)]);

end