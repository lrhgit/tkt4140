function yd = fcndeflect(x,y)
% Used by program Cantilev1 - and 2
global alpha2 
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = -alpha2*cos(y(1));
yd(3) = sin(y(1));