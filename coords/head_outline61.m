function [] = head_outline61(plot_electrodes)
% Create a circle using rectangle command
rectangle('Position',[-.95,-.95,1.9,1.9],'Curvature',[1,1],'LineWidth',2)
axis([-1.05,1.05,-1.00,1.05]);
axis square
hold on

% Create a 'nose' using two lines joining at center ('line' command
line([-.12; 0; .12], [.945; 1.05; .945], 'color', 'k','LineWidth',2);   % Nose

% Create Right-side and Left-side Ears
Xr = [0.9450 0.9730 0.9810 0.9929 1.0049 1.0030 1.0100 0.9930 0.9730 0.9300];
Xl = Xr .* -1.0;
Y  = [0.0555 0.0775 0.0783 0.0746 0.0555 -0.0055 -0.0932 -0.1313 -0.1384 -0.1199];
Yc = Y .* 1.5;
line(Xr,Yc,'color','k','LineWidth',2)
line(Xl,Yc,'color','k','LineWidth',2)
%     colorbar('vert');

% Create dots for each electrode position
if plot_electrodes
    load('/export/home/mike/matlab/origin/coords/elec_positions.mat');
    plot(elec_pos(:, 1),elec_pos(:, 2),'k.','MarkerSize',10)
end