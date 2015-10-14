%load location file
%coords_filepath='A:\matlab\coords\31chans_ns.mat';
seed_chanlocs=load('/export/home/mike/matlab/origin/coords/31chans_ns.mat');
seed_chanlocs=getfield(seed_chanlocs,'chan_locs');
n_channels=length(seed_chanlocs);

%specify electrodes to pair with all others
% FZ, CZ, PZ
electrode_seeds=[17 18];
n_seeds=length(electrode_seeds);

%generate pairs
pairs=zeros(n_seeds*(n_channels-1),2);
for seed=1:length(electrode_seeds)
    pairs( ((n_channels-1)*(seed-1)+1):((n_channels-1)*seed), : ) = ...
        [ ones(n_channels-1,1)*electrode_seeds(seed) , setdiff(1:n_channels,electrode_seeds(seed))' ];
end

clear coords_filepath seed_chanlocs n_channels electrode_seeds n_seeds seed