function plotlabel(xlabel,ylabel)

font_size=10;

%"xlabel"
ax = axes('Units', 'Normal', 'Position', [.075 .0001 .85 .0001], ...
		'Visible', 'off');
set(get(ax, 'Title'), 'Visible', 'on')
title(xlabel, 'FontSize', font_size);    

%"ylabel"
ax = axes('Units', 'Normal', 'Position', [.062 .05 .025 .5], ...
		'Visible', 'off');
set(get(ax, 'Title'), 'Visible', 'on')
yt=title(ylabel, 'FontSize', font_size);    
set(yt,'rotation',90);

end