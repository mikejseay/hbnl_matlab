time = 48;
freq = 12;
cond = 1;
pair = 54;

n = meanx(n_trials_all,[]);
ispc_val = meanx(cohdata(time,freq,cond,pair,s_inds_g(:,1)),[]);
phi_d = angle( meanx( exp(1i*phidata(time,freq,cond,pair,s_inds_g(:,1))) , [] ) );

gv = n .* ispc_val .* exp( (-1.*phi_d.^2) ./ (4.*pi./n) ) .* ( (2./n) .^ 0.5 );