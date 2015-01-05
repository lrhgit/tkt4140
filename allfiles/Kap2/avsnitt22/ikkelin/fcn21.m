function dydx = fcn21(x,y)
global a;
dydx = zeros(2,1);
dydx(1) = y(2);
dydx(2) = a*y(1)^2;
