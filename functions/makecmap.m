function cmap=makecmap(lims, center, colors, n_shades)
% makes a color map with white aligned with center and colorbrewer gradients
% toward positive (red) and negative limits (blue) indicated by 2-element
% vector lims. if the center is not in the middle of lims,
% uses a gradient from lims(1) (white) to high (green).

% example: cmap = makecmap([-.2 .2], 0, 'redblue', 64);

if nargin < 4
    n_shades = 256;
end
if nargin < 3
    colors = [];
end
if nargin < 2
    center = 0;
end

% color override
if ~isempty(colors)
    switch colors
        case 'redblue'
            poscolor='Reds';
            negcolor='Blues';
        case 'purpleorange'
            poscolor='Purples';
            negcolor='Oranges';
    end
else
    if center==0
        poscolor = 'Reds';
        negcolor = 'Blues';
    elseif center == 1
        poscolor = 'Purples';
        negcolor = 'Oranges';
    end 
end

if lims(1) < center && lims(2) > center %(polar data)
    
    [~, white_ind]=min(abs( linspace(lims(1),lims(2),n_shades) - center ));
    if (n_shades - white_ind) > white_ind
        cmap_pos = cbrewer('seq',poscolor,n_shades-white_ind);
        cmap_neg = flipud(cbrewer('seq',negcolor,n_shades-white_ind));
        cmap_neg(1 : (n_shades - 2 * white_ind + 1), :)=[];
        %cmap_pos(end,:)=[];
        cmap=[cmap_neg; [1 1 1]; cmap_pos];
    else
        cmap_pos=cbrewer('seq', poscolor, white_ind);
        cmap_neg=flipud(cbrewer('seq', negcolor, white_ind));
        cmap_pos((n_shades+1-white_ind) : end, :)=[];
        cmap_neg(end, :)=[];
        cmap=[cmap_neg; [1 1 1]; cmap_pos];
    end
    
else
       
    cmap=cbrewer('seq', 'Greens', n_shades);
    
end
        



end