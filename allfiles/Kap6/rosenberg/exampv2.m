%program exampv2
% Motstrømsvarmeveksler
% Eksempel 3-1 i Rosenberg, side 44
%
clear; close
clearvars -global b1 b4 b5 b6 c4 c6;
global b1 b4 b5 b6 c4 c6;
nx = 20; % antall xskritt
nt = 20; % antall tidskritt
dx = 1/nx;
nxm1 = nx - 1;
nxp1 = nx + 1;
u = zeros(nx,1); d1 = u; d2 = u;
v = zeros(nx,1);
cc1 = 0.05;
cc2 = cc1;
bb2 = 0.5;
b1 = 1 + 4*nx/cc1;
b5 = 4*nx/cc1 - 1;
b4 = -(1 + 2*nx*(bb2 + 1)/cc2);
c6 = 1 - 2*nx*(bb2 + 1)/cc2;
c4 = 2*nx*(bb2 - 1)/cc2 - 1;
b6 = 2*nx*(bb2 - 1)/cc2 + 1;
dt = dx;
x = (0:dx:1)';
t = 0;
% Utskrift
for n = 1: nt
    t = t + dt;
    d1(1) = -1 + b5 -u(1) + v(2) + v(1);
    d2(1) = - 2 - u(1) + b6*v(1) + c6*v(2);
    for i = 2:nxm1
        d1(i) = b5*u(i-1) - u(i) + v(i+1) + v(i);
        d2(i) = -u(i) - u(i-1) + b6*v(i) + c6*v(i+1);
    end
     d1(nx) = b5*u(nxm1) - u(nx) + v(nx);
     d2(nx) = -u(nx) - u(nxm1) + b6*v(nx);
    [u,v] = bitris2(d1,d2);
end
u =[1;u];
v = [v;0];
% plot(x,u,x,v);
%     


