% quick test to try a motif for plotting separate relational pair hypotheses

hyps=unique(hyp_inds);
hyp_spdims=numSubplots(length(hyps));

figure;
subplot_dummy=0;
for hyp=hyps'
    subplot_dummy=subplot_dummy+1;
    subplot(hyp_spdims(1),hyp_spdims(2),subplot_dummy);
    inds2p=find(hyp_inds==hyp);
    plot_pairs(pairs(inds2p,:),hyp_inds(inds2p),chan_locs);
end