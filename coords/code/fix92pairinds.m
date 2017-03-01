opt.coherence_inds=zeros(92,1);
opt.coherence_inds(1:13) =1;  %frontal
opt.coherence_inds(14:26)=2; %parietal
opt.coherence_inds(27:37)=3; %occipital
opt.coherence_inds(38:59)=4; %frontal-parietal
opt.coherence_inds(60:77)=5; %frontal-occipital
opt.coherence_inds(78:92)=6; %parieto-occipital

opt.coherence_indlbls={'frontal','parietal','occipital','frontal-parietal', ...
    'frontal-occipital','parietal-occipital'};

%%

hyp_inds=zeros(92,1);
hyp_inds(1:13) =1;  %frontal
hyp_inds(14:26)=2; %parietal
hyp_inds(27:37)=3; %occipital
hyp_inds(38:59)=4; %frontal-parietal
hyp_inds(60:77)=5; %frontal-occipital
hyp_inds(78:92)=6; %parieto-occipital

hyp_indlbls={'frontal','parietal','occipital','frontal-parietal', ...
    'frontal-occipital','parietal-occipital'};