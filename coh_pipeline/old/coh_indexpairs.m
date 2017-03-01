function scl = coh_indexpairs(opt, scl)

% create pair indexes based on 1.) seeds 2.) regions 3.) inter vs. intra

n_pairs = length(opt.coherence_pairs);
