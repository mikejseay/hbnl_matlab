function set_print_size(width, height, linewidth, orientation)
	set(gcf, 'PaperPositionMode', 'manual');
	set(gcf, 'PaperUnits', 'centimeters');
	if nargin > 3 && ~isempty(orientation)
		set(gcf, 'PaperOrientation', orientation);
	end	
	set(gcf, 'PaperPosition', [1, 1, width, height]);
	if nargin > 2 && ~isempty(linewidth)
		lineobj = findobj('type', 'line');
		set(lineobj, 'linewidth', linewidth);
	end


