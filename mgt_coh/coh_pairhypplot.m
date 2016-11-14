function coh_pairhypplot(pairs,chan_locs,hyp_inds,hyp_indlbls,clickable,direction_spec)
% plot separate relational pair hypotheses

if nargin<5
    clickable=false;
end
if nargin<6
    direction_spec=1;
end

hyps=unique(hyp_inds);
hyp_spdims=numSubplots(length(hyps));

figure;
subplot_dummy=0;
for hyp=hyps'
    subplot_dummy=subplot_dummy+1;
    subplot(hyp_spdims(1),hyp_spdims(2),subplot_dummy);
    inds2p=find(hyp_inds==hyp);
    plot_pairs(pairs(inds2p,:),chan_locs,hyp_inds(inds2p),clickable,direction_spec);
    title(hyp_indlbls{hyp});
end

end