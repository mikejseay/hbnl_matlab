% build a channel locations file from scratch

load('61chans.mat')
chan_locs=readlocs( '/export/home/mike/matlab/eeglab13_4_4b/sample_locs/Standard-10-10-Cap47.ced');

for chan=1:61
    %insert the labels and XYZ coords
    chan_locs(chan).labels=chan_names{chan};
    chan_locs(chan).Y=-chan_coords(chan,1);
    chan_locs(chan).X=chan_coords(chan,2);
    chan_locs(chan).Z=chan_coords(chan,3);    
end

%get rid of everything else
%fields_removed={'theta','radius','sph_theta','sph_phi','sph_radius','sph_theta_besa','sph_phi_besa','type'};
%for field=1:length(fields_removed)
%    chan_locs=rmfield(chan_locs,fields_removed{field});
%end
