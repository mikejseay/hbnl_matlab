function colorscale_plot(lims, cmap, mid_pos, n_colors, n_ticks, precision, orientation)
% wrapper for colorscale to plot onto separate figure

if nargin < 7
    orientation = 'vert';
end
if nargin < 6
    precision = 2;
end
if nargin < 5
    n_ticks = 5;
end
if nargin < 4
    if nargin > 1
        n_colors = size(cmap, 1);
    else
        n_colors = 256;
    end
end
if nargin < 3
    mid_pos = 0.5;
end
if nargin < 2
    cmap = parula(n_colors);
end

switch orientation
    case 'vert'
        cb_pos=[mid_pos 0.1 0.05 0.8];
    case 'horiz'
        cb_pos=[0.1 mid_pos 0.8 0.05];
end

h=colorscale([1 n_colors], lims, range(lims)/n_ticks, orientation, precision, 'Position', cb_pos);
set(h,'FontSize',20);

colormap(h,cmap);

end