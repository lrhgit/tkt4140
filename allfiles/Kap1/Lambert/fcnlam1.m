function dydx = fcnlam1(x,y)
dydx = zeros(3,1);
dydx(1) = - 21*y(1) + 19*y(2) - 20*y(3);
dydx(2) = 19*y(1) - 21*y(2) + 20*y(3);
dydx(3) = 40*y(1) - 40*y(2) - 40*y(3);

