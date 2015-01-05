function yd = beam(x,y)
% Used by program cantlevplot
global alpha2
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = -alpha2*cos(y(1));
