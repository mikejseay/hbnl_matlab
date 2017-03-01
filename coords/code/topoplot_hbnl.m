function topoplot_hbnl(data2plot, subtitle, scale_min, scale_max, contour_levels, ...
    outline_width, dots_size, ColorBar_YesHorz1_YesVert2_No0,cmap,plot_electrodes)
% Topo_HeadPlot_MATvar_Single_Plain(data2plot, subtitle, scale_min, scale_max,
% contour_levels, outline_width, dots_size, ColorBar_YesHorz1_YesVert2_No0)
% [This program plots the scalp topography (i.e. HeadPlot) for the given data]
% {***Program compiled by Chella Kamarajan***)

% data2plot -> data loaded as matlab variable
% scale_min -> minimum value for the color scale
% scale_max -> maximum value for the color scale
% contour_levels -> number of contour levels (color grades)
% outline_width -> line thickness for the head outlines (circles, ears, and nose)
% dots_size -> size of elctrode markers (dots) in font sizes
% ColorBar_YesHorz1_YesVert2_No0 -> Color bar along with headplot: Yes-Horizontal '1', Yes-Vertical '2', No '0'

%% Data and messages

if size(data2plot,1) >= 61 && size(data2plot,2) > 1
    fprintf('Data has multiple columns...!! Program will use the first data column')
end

if size(data2plot,1) < 61
    error('Data column has less than 61 values...!! Program needs 61 values...!!!')
elseif size(data2plot,1)==64
    fprintf('Data column has 64 values...!! Program will use only 61 values of scalp sites (out of 64)...!!!')
    data2plot = data2plot([1:31,33:62],:);
elseif size(data2plot,1)>64
    error('Data column has more than 64 values...!! Program needs 61 or 64 values...!!!')
end

%% Input Arguments

if nargin < 1
    error('Data has not been supplied or identified...!!! -> Please check the input data again!!!')
end

if nargin < 2
    subtitle = [];
end

if nargin < 3
    scale_min = min(min(data2plot));
end

if nargin < 4
    scale_max = max(max(data2plot));
end

if nargin < 5
    contour_levels = 15;
end

if nargin < 6
    outline_width = 2;
end

if nargin < 7
    dots_size = 5;
end

if nargin < 8
    ColorBar_YesHorz1_YesVert2_No0 = 1;
end

if nargin < 9
    cmap=jet(64);
end

if isempty(scale_min)
    scale_min = min(min(data2plot));
end

if isempty(scale_max)
    scale_max = max(max(data2plot));
end


%% Making a 2D grid
xy_scale = -1.05:0.025:1.05;
[xi, yi] = meshgrid(xy_scale);
z_grid = xi + sqrt(-1) * yi;
abs_z_grid = z_grid .* conj(z_grid);
z_template = find(abs_z_grid <= 1.0);
xi_circ = xi(z_template);
yi_circ = yi(z_template);
z_mat = zeros(size(abs_z_grid));
z_mat = z_mat + NaN;

%% Create Coordinates Data
load('/export/home/mike/matlab/origin/coords/electrode_coordinates_61.mat');

%% Creating Head Plots

Z = griddata(coord61(:, 1), coord61(:, 2), data2plot, xi_circ, yi_circ, 'cubic'); %#ok<*GRIDD>
z_mat(z_template) = Z;

if size(data2plot,1)==61
    head_outline61(plot_electrodes);
end

[~, H] = contourf(xi, yi, z_mat, contour_levels);
set(H, 'linestyle', 'none');
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
hold on;                % Further feature(s) on this to be added

% Create a circle using rectangle command
rectangle('Position',[-.95,-.95,1.9,1.9],'Curvature',[1,1],'LineWidth',outline_width)
axis([-1.05 1.05 -0.95 1.05]);
axis square
axis off
hold on

% Create a 'nose' using two lines joining at center ('line' command)
line([-.12; 0; .12], [.945; 1.05; .945], 'color', 'k','LineWidth',outline_width);   % Nose

% Create right-side and left-side ears
Xr = [0.9450 0.9730 0.9810 0.9929 1.0049 1.0030 1.0100 0.9930 0.9730 0.9300];
Xl = Xr .* -1.0;
Y  = [0.0555 0.0775 0.0783 0.0746 0.0555 -0.0055 -0.0932 -0.1313 -0.1384 -0.1199];
Yc = Y .* 1.5;
line(Xr,Yc,'color','k','LineWidth',outline_width)
line(Xl,Yc,'color','k','LineWidth',outline_width)

% Create dots for each electrode position
if plot_electrodes
    elec_pos = load('Elec_Positions.csv');
    plot(elec_pos(:,1),elec_pos(:,2),'k.','MarkerSize',dots_size)
end
    
caxis([scale_min scale_max]); axis off;
colormap(cmap);

if ColorBar_YesHorz1_YesVert2_No0 == 1
    colorbar('SouthOutside','FontSize',7);
elseif ColorBar_YesHorz1_YesVert2_No0 == 2
    colorbar('EastOutside','FontSize',7);
elseif ColorBar_YesHorz1_YesVert2_No0 == 0
    colorbar('off');
end

if ~isempty(subtitle)
    title(subtitle,'FontSize',12,'FontWeight','normal','Interpreter','none');
end
end

%% END