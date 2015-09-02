%     Even though John and Duane have given you solutions to your problem,
%     Daniel, I cannot resist giving one in a more condensed form. If (x1,y1)
%     and (x2,y2) are the beginning and end points, respectively, and r the
%     radius of the desired arc which is to bend in a clockwise (default)
%     direction. % Original: plot(x,y,'y.',x1,y1,'r*',x2,y2,'b*')
%
%     ---by   Roger Stafford


function Draw_Arc_Clockwise(Start_Point, End_Point, ColorTriplet, SynchLineWidth, Direction)
% Provide Start_Point and End_Point as (x1,y1) and (x2,y2); 
% e.g., Draw_Arc_Clockwise([0,0.40], [0,-0.4], 'b-', 'bo', 'bo', 3);

x1 = Start_Point(1);
y1 = Start_Point(2);
%z1 = Start_Point(3);
x2 = End_Point(1);
y2 = End_Point(2);
%z2 = End_Point(3);


d = sqrt((x2-x1)^2+(y2-y1)^2); % Distance between points
if Direction==1
    a = atan2(x2-x1,-(y2-y1)); % Perpendicular bisector angle (for clockwise arc)
else
    a = atan2(-(x2-x1),y2-y1); % Perpendicular bisector angle (for anticlockwise arc)
end
r = d;                       % Radius of the circle (length of the circle as default)
b = asin(d/2/r); % Half arc angle
c = linspace(a-b,a+b); % Arc angle range
e = sqrt(r^2-d^2/4); % Distance, center to midpoint
x = (x1+x2)/2-e*cos(a)+r*cos(c); % Cartesian coords. of arc
y = (y1+y2)/2-e*sin(a)+r*sin(c);

Fig_Handle = plot(x,y,'Color',ColorTriplet);
set(Fig_Handle, 'LineWidth',SynchLineWidth); 
% axis equal

end



