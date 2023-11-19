function addMrkr(x,y,R,cL,cR)
% x,y are row vectors of line points
% R is radius
% cL,cR are colors of left/right semicircles
% left
th = linspace(pi/2, 3*pi/2)';
xL = R*cos(th) + x;
yL = R*sin(th) + y;
hold on
patch(xL,yL,cL,'LineStyle','none');

% % right
th = linspace(pi/2, -pi/2)';
xR = R*cos(th) + x;
yR = R*sin(th) + y;
h = patch(xR,yR,cR,'LineStyle','none');
hold off
end