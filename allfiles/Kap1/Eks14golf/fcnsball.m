function dydt = fcnsball(t,y)
% Used in program sballrk4c
% Called by RK4C
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
