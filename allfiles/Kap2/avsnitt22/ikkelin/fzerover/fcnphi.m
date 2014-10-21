function phi = fcnphi(s)
% Brukes av programmet Eks22
xspan = [0.0 1.0];
y0 = [4.0 s];
options = odeset('RelTol',1.0e-5);
[x,y] = ode45(@fcn22,xspan,y0,options);
phi = y(end,1) - 1;   