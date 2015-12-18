function colorscale_plot(lims, cmap, horiz_pos, n_colors, n_ticks, precision)
% wrapper for colorscale to plot onto separate figure

if nargin < 6
    precision = 2;
end
if nargin < 5
    n_ticks = 5;
end
if nargin < 4
    n_colors = 256;
end
if nargin < 3
    horiz_pos = 0.5;
end
if nargin < 2
    cmap = parula(n_colors);
end

cb_pos=[horiz_pos 0.1 0.05 0.8];

h=colorscale([1 n_colors], lims, range(lims)/n_ticks, 'vert', precision, 'Position', cb_pos);
set(h,'FontSize',20);

colormap(h,cmap);

end