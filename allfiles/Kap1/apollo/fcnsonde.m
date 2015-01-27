function dydt = fcnsonde(t,y)
dydt = zeros(size(y));
global mu lam;
r1 = sqrt((y(1) + mu)^2 + y(2)^2);
r2 = sqrt((y(1) - lam)^2 + y(2)^2);
r13 = r1^3; r23 = r2^3;
dydt(1) = y(3);
dydt(2) = y(4);
dydt(3) = y(1) + 2*y(4) - lam*(y(1) + mu)/r13 - mu*(y(1) - lam)/r23;
dydt(4) = y(2) - 2*y(3) -(lam/r13 + mu/r23)*y(2);
