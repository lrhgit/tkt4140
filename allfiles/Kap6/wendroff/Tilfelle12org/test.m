%program test
% Motstrømsvarmeveksler
% Bygger på eksempel 3-1 i Rosenberg, side 44
%
clear; close
clearvars -global b1 b4 c4 ;
global b1 b4 c4;
nx = 50; % antall xskritt
nt = 60; % antall tidskritt
N = 5; % Antall skritt for økning av u0
dx = 1/nx;
nxm1 = nx - 1;
nxp1 = nx + 1;
u = zeros(nx,1); d1 = u; d2 = u;
v = zeros(nx,1);
Tyr = 100; Tir = 30;
Tdiff = Tyr - Tir;
alfy = 0.9;
alfi = 0.2;
b2 = 1.0; % b2 = wi/wy;
b1 = 1 + 4*nx/alfy;
b5 = 4*nx/alfy - 1;
b4 = -(1 + 2*nx*(b2 + 1)/alfi);
c6 = 1 - 2*nx*(b2 + 1)/alfi;
c4 = 2*nx*(b2 - 1)/alfi - 1;
b6 = 2*nx*(b2 - 1)/alfi + 1;
dt = dx;
x = (0:dx:1)';
t = 0;
% Utskrift
for n = 1: nt
    t = t + dt;
    if n <= N
        f1 = n/N;
        f2 = (n+1)/N;
        % Lineær økning av u0
        u0n = f1;
        u0np1 = f2;
        % Økning etter 2.gradspolynom
%         u0n = f1*(2 - f1);
%         u0np1 = f2*(2 - f2);
        % Økning etter 3.gradspolynom
%         u0n = f1*(3 - f1^2)/2;
%         u0np1 = f2*(3 - f2^2)/2;      
    else
        u0n = 1;
        u0np1 = 1;
    end
    d1(1) = -u0np1 + b5*u0n - u(1) + v(2) + v(1);
    d2(1) = - (u0n + u0np1) - u(1) + b6*v(1) + c6*v(2);
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
Ty = Tdiff*u + Tir;
Ti = Tdiff*v + Tir;
plot(x,Ty);
%     


