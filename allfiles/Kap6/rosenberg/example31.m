%program example31
% Motstrømsvarmeveksler
% Eksempel 3-1 i Rosenberg, side 44
%
clear; close
nx = 20; % antall xskritt
nt = 20; % antall tidskritt
dx = 1/nx;
nxm1 = nx - 1;
nxp1 = nx + 1;
u = zeros(nx,1); g1 = u; g2 = u;
el2 = u; el4 = u;
v = zeros(nxp1,1);
cc1 = 0.05;
cc2 = cc1;
bb2 = 0.5;
t1 = 4*nx/cc1;
b1 = 1 + t1;
b5 = t1 - 1;
t2 = 2*nx*(bb2 + 1)/cc2;
b4 = -(1 + t2);
c6 = 1 - t2;
t3 = 2*nx*(bb2 - 1)/cc2;
c4 = t3 - 1;
b6 = t3 + 1;
dt = dx;
x = (0:dx:1)';
t = 0;
u0 = 1;
elt = b1*c4 + 1;
% Utskrift
for n = 1: nt
    t = t + dt;
    em = b1*b4 - 1;
    el2(1) = (c4 - b4)/em;
    el4(1) = elt/em;
    d1 = -1 + b5 -u(1) + v(2) + v(1);
    d2 = - 2 - u(1) + b6*v(1) + c6*v(2);
    g1(1) = (b4*d1 + d2)/em;
    g2(1) = (b1*d2 - d1)/em;
    for i = 2:nx
        bt2 = - 1 - el2(i-1);
        bt4 = b4 - el2(i-1);
        em = b1*bt4 - bt2;
        el2(i) = (-bt4 - c4*bt2)/em;
        el4(i) = elt/em;
        d1 = b5*u(i-1) - u(i) + v(i+1) + v(i);
        d2 = -u(i) - u(i-1) + b6*v(i) + c6*v(i+1);
        g1(i) = (bt4*(d1 - g1(i-1)) - bt2*(d2 - g1(i-1)))/em;
        g2(i) = (b1*(d2 - g1(i-1)) - d1 + g1(i-1))/em;
    end
    u(nx) = g1(nx);
    v(nx) = g2(nx);
    for j = 1: nxm1
        i = nx - j;
        u(i) = g1(i) - el2(i)*v(i+1);
        v(i) = g2(i) - el4(i)*v(i+1);
    end
end
u =[1;u];
plot(x,u,x,v);
    


