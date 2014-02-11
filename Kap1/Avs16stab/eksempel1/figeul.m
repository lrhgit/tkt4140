% Program figeul
% Løser eksempel 1 i avsnitt 1.6.2
% av Eulers metode for dt = 0.011
% for å demonstrere instabilitet.
% Løsningsintervall : 0 <= t <= 0.22
% Plotter løsning og sammenligner med den analytiske

% Ligning : x''(t)= l1*x'(t) + l*x(t)
%           x(0) = 1, x'(0) = 0 = v(0)
%           l1 = -2*c, l = -1
%           Bruker c = 100 -> l1 = -200
clear; clf;
c = 100;
l1 = -2*c; l = -1;
tmax =  0.22;
dt = 0.011;
nmax = round(tmax/dt) + 1; % Number of time-steps
% Initializing of av vectors
t = zeros(nmax,1); x = t; v = t;

% Initial values
x(1) = 1; v(1) = 0;  t(1) = 0;

% === Euler's method ===
for n = 1 : nmax - 1
    x(n+1) = x(n) + dt*v(n);
    v(n+1) = (1.0 + l1*dt)*v(n) + dt*l*x(n);
    t(n+1) = dt*n;
end 

% === Analytisk løsning i endepunktet
d = sqrt(c^2 - 1);
b1 = -c + d; b2 = -c -d; b12 = 2*d;
A = -b2/b12; B = b1/b12;
xa = A*exp(b1*t) + B*exp(b2*t) ;
plot(t,x,'k',t,xa,'k','LineWidth',1)
shg

grid on
FS = 'FontSize'; FW = 'FontWeight';
st = sprintf('Eksempel 1 med \\Deltat = 0.009 og 0.011');
title(st,FS,14,FW,'Bold');
xlabel('t',FS,14,FW,'Bold');
ylabel('x',FS,14,FW,'Bold');