function yd = fcn22(x,y)
%--- Brukes av programmet eks22 og plot22
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = 3*y(1)^2/2;