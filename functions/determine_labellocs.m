function [locs]=determine_labellocs(dim)

v=axis;
xmin=v(1);
xmax=v(2);
ymin=v(3);
ymax=v(4);
width=xmax-xmin;
height=ymax-ymin;

if strcmp(dim,'row')
    locs= [ xmin - (0.7*width), ymax - (0.4*height) ];
elseif strcmp(dim,'column')
    locs= [ xmin + (0.1*width), ymax - (1.15*height) ];
end


end