function [pairs,seed_inds,seed_indlbls]=pair_with_all_others(seed_chanlocs, seed_electrodes)
% pairs an array of electrodes (seed_electrodes) with all other electrodes
% in the seed_chanlocs channel locations structure, generating a "seed
% coherence" hypothesis of pairs

% load location file
if ~isstruct(seed_chanlocs)
    seed_chanlocs=load(seed_chanlocs);
    seed_chanlocs=getfield(seed_chanlocs,'chan_locs');
end

n_channels=length(seed_chanlocs);
chan_labels={seed_chanlocs.labels};

% if seed_electrodes are cell array of strings (channel names, interpret)
if ~all(ismember(seed_electrodes,chan_labels))
    error('Channel name incorrectly specified, use capital letters')
else
    if iscell(seed_electrodes)
        chanind=zeros(length(seed_electrodes),1);
        for seedind=1:length(seed_electrodes)
            chanind(seedind)=find(strcmpi(chan_labels,seed_electrodes{seedind}));
        end
        seed_electrodes=chanind;
    elseif ischar(seed_electrodes)
        seed_electrodes=find(strcmpi(chan_labels,seed_electrodes));
    end
end

%specify electrodes to pair with all others
n_seeds=length(seed_electrodes);

%generate pairs & indices
pairs=zeros(n_seeds*(n_channels-1),2);
seed_indlbls=cell(n_seeds,1);
for seed=1:n_seeds
    pairs( ((n_channels-1)*(seed-1)+1):((n_channels-1)*seed), : ) = ...
        [ ones(n_channels-1,1)*seed_electrodes(seed) , setdiff(1:n_channels,seed_electrodes(seed))' ];
    seed_inds( ((n_channels-1)*(seed-1)+1):((n_channels-1)*seed), : ) = ...
        ones(n_channels-1,1)*seed;
    seed_indlbls{seed}=seed_chanlocs(seed_electrodes(seed)).labels;
end

end