function dydx = fcntank1(x,y) 
% Used by program tank1
global beta4;
dydx = zeros(size(y));
dydx(1) = y(2);
dydx(2) = y(3);
dydx(3) = y(4);
dydx(4) = -4*beta4*(y(1) + 1 - x);
