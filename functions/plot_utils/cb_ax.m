function cb_ax(dist, width, cb_lim, n_colors, ticks, orientation)

cb_pos = get(gca,'Position') + [dist 0 0 0];
cb_pos(3) = width; %set width
colorscale([1 n_colors], cb_lim, range(cb_lim)/ticks, orientation, ...
    'Position',cb_pos);
end