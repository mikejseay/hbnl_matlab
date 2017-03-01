function maskersp=coh_checksig(P,t_start_b,t_end_b,alpha)

% P = meanx(wave_totdata(:,7,:,2,:),[1 3])';
% [~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
% [~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));
% alpha=.05;

formula='arg1';
n_perms=200;

fdr_bool=true;

[ Pboot Pboottrialstmp Pboottrials] = bootstat(P, formula, 'boottype', 'shuffle', ...
    'label', 'ERSP', 'bootside', 'both', 'naccu', n_perms, ...
    'basevect', t_start_b:t_end_b, 'alpha', alpha, 'dimaccu', 2 );
    
%bootside should be 'upper' for ITC

Pboottrials=Pboottrials';
exactp_ersp = compute_pvals(P, Pboottrials);
alphafdr = fdr(exactp_ersp, alpha);

if fdr_bool
    maskersp = exactp_ersp <= alphafdr;
else
    maskersp = exactp_ersp <= alpha;
end

end