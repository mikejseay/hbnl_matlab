function inds = balance_dists1(a, b)

% given 2 vectors (a and b) of matched histogram bin counts, returns indices
% to remove from the original sorted vector of b which will best balance
% the distributions i.e. b has more observations total and the
% distributions don't match

n_bins = length(a);

a_cum = cumsum(a);
b_cum = cumsum(b);

diff_tot = sum( subplus( b - a ) );
diff_tot2 = b_cum(end) - diff_tot;
diff_tot3 = a_cum(end) - diff_tot2;
numbins2fix = sum( subplus( b - a ) > 0 );
lessen_factor = floor(diff_tot3 ./ numbins2fix) + 3; % conservative

inds = [];
for bin = 1:n_bins
    
    if b(bin) > a(bin)
        
        amt2remove = b(bin) - a(bin) - lessen_factor;
        
        if amt2remove > 0
            inds = [ inds, sort(randperm( b(bin), amt2remove )) + b_cum(bin - 1) ];
        end
        
    end
    
end