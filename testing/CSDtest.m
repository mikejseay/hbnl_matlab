
% specify the channel locations in eeglab format
chanlocs_file='/export/home/mike/matlab/origin/coords/61chans_ns.ced';

% convert locations from .ced to .csd for CSD
ConvertLocations(chanlocs_file);
csd_file=[chanlocs_file(1:end-3),'csd'];
mat_file=[chanlocs_file(1:end-3),'mat'];

% load the channel locations into memory and create an M matrix for CSD
load(mat_file);
E={chan_locs.labels}';
M=ExtractMontage(csd_file,E);

% make sure the montage looks right
MapMontage(M);

% calculate the G and H matrices using default spline rigidity
[G,H]=GetGH(M);