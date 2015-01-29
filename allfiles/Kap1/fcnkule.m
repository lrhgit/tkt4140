function dy = fcnkule(x,y)
% Brukes av rkule
n = length(y);
dy = zeros(n,1);
global A B C d nu Re;
v = y(2);
va = abs(v);
Re = va*d/nu;
CD = CDkule(Re);
dy(1) = v;
dy(2) = (B - C*v*va*CD)/A;