function dydt = fcngball(t,y)
% Used in program golfrk4c
% Called by RK4C
dydt = zeros(size(y));
global g C d nu Re vfx vfy nrpm;
vrx = y(3) - vfx;
vry = y(4) - vfy;
vr = sqrt(vrx^2 + vry^2);
Re = vr*d/nu;
[cd , cl] = cdcldata(vr,nrpm);
f1 = C*vr*cd; f2 = C*vr*cl;
dydt(1) = y(3);
dydt(2) = y(4);
dydt(3) = - f1*vrx - f2*vry;
dydt(4) = - g + f2*vrx -f1*vry;
