function dydt = kule2d45(t,y)
% Used by the program golf45
% Called by ode45
% Dragdata from function CDkule
dydt = zeros(size(y));
global g C d nu Re vfx vfy;
vrx = y(3) - vfx;
vry = y(4) - vfy;
vr = sqrt(vrx^2 + vry^2);
Re = vr*d/nu;
CD = CDkule(Re);
f = C*vr*CD;
dydt(1) = y(3);
dydt(2) = y(4);
dydt(3) = - f*vrx;
dydt(4) = - g -f*vry;
