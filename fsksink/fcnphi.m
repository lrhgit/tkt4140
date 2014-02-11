function phi = fcnphi(s)
% Called from program SINK
global x1; %ksi-inf
xspan = [0 x1];
g0 = 2*(sqrt(3) - 1);
y0 = [g0 ; -4/3 ; s];
tol = 1.0e-9;
options = odeset('RelTol',tol,'AbsTol',[tol tol tol*1.0e-4]);
[x,y] = ode15s(@fcnsink,xspan,y0,options);
phi = y(end,2) -1;

