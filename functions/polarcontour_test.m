%test
% Create polar data
[r,t] = meshgrid(0:.1:5,0:pi/30:(2*pi));
z = r - t;
% Convert to Cartesian
x = r.*cos(t);
y = r.*sin(t);
h = polar(x,y);
hold on;
contourf(x,y,z);
% Hide the POLAR function data and leave annotations
set(h,'Visible','off')
% Turn off axes and set square aspect ratio
axis off
axis image