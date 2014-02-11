function dy = kule2d(x,y)
% Brukes i programmet Shoot
% Kalles av RK4C
dy = zeros(size(y));
global g C d nu Re vfx vfy;
vrx = y(3) - vfx;
vry = y(4) - vfy;
vr = sqrt(vrx^2 + vry^2);
Re = vr*d/nu;
CD = CDkule(Re);
f = C*vr*CD;
dy(1) = y(3);
dy(2) = y(4);
dy(3) = - f*vrx;
dy(4) = - g -f*vry;
