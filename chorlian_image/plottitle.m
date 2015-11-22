function [ax, h] = plottitle(text, pos, font_size)
%
%Centers a title over a group of subplots.
%Returns a handle to the title and the handle to an axis.
% [ax, h]= subtitle(text)
% returns handles to both the axis and the title.
% ax= subtitle(text)
% returns a handle to the axis only.
if nargin == 1 || isempty(pos)
	pos = 1;
end
if nargin < 3 || isempty(font_size)
	font_size = 12;
end
if pos == -1
	ax = axes('Units', 'Normal', 'Position', [.075 .005 .85 .05], ...
		'Visible', 'off');
elseif pos == 1
	% original axis vector = [.075 .92 .85 .05],
	ax = axes('Units', 'Normal', 'Position', [.075 .835 .82 .10], ...
			'Visible', 'off');
% 	if font_size < 10
% 		ax = axes('Units', 'Normal', 'Position', [.075 .85 .85 .12], ...
% 			'Visible', 'off');
% 	else
% 		ax = axes('Units', 'Normal', 'Position', [.075 .80 .85 .15], ...
% 			'Visible', 'off');		
% 	end
else
	ax = []; h = []; return;
end
text = strrep(text, '_', ' ');
set(get(ax, 'Title'), 'Visible', 'on')
title(text, 'FontSize', font_size);
if (nargout < 2)
return
end
h= get(ax, 'Title');
