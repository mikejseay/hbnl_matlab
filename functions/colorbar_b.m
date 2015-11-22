function [ax, h] = colorbar_b()
% color bar but easier to size
ax = axes('Units', 'Normal', 'Position', [.475 .25 .5 .5], ...
    'Visible', 'off');
colorbar;
h= get(ax, 'Title');