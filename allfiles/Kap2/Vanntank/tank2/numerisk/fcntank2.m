function dydx = fcntank2(x,y) 
% Used by program tank2
global beta4 alfa;
z = 1 - alfa*x;
dydx = zeros(size(y));
dydx(1) = y(2);
dydx(2) = y(3);
dydx(3) = y(4);
temp = (6*alfa/z)*y(4) -(6*alfa^2/z^2)*y(3);
dydx(4) = temp -4*beta4*y(1)/z^2 - 4*beta4*(1 - x)/z^3;
