% Program figeul
% Løser eksempel 2 i avsnitt 1.6.2
% av Eulers metode for dt = 1.55e-5
% for å demonstrere instabilitet.
% Løsningsintervall : 0 <= t <= 0.22
% Plotter løsning og sammenligner med den analytiske

% Ligning : y''(t)= l1*y'(t) + l*y(t) + c
%           y(0) = 1, y'(0) = 0 = v(0)
%           l1 = -129600, l = 98696, c = 9869.6
clear; clf;
l1 = -129600; l = 98696 ; c = 9869.6 ;
tmax =  1.2e-3;
dt = 1.58e-5;
nmax = round(tmax/dt) + 1; % Number of time-steps
% Initializing of av vectors
t = zeros(nmax,1); y = t; y = t;

% Initial values
y(1) = 0; v(1) = 0;  t(1) = 0;

% === Euler's method ===
for n = 1 : nmax - 1
    y(n+1) = y(n) + dt*v(n);
    v(n+1) = (1.0 + l1*dt)*v(n) + dt*(c + l*y(n));
    t(n+1) = dt*n;
end 

% === Analytisk løsning 

% b1 = 0.761538735021; b2 = -129600.7615387;
% A = 0.09999941239986; B = 5.87600143067797e-7;
d = sqrt(l1^2 +4*l);
b1 = (l1 + d)*0.5; b2 = (l1 - d)*0.5 ; b12 = -b1/b2;
A = 0.1/(1 + b12); B = b12*A;
ya = A*exp(b1*t) + B*exp(b2*t) - 0.1;
plot(t,y,'k',t,ya,'k','LineWidth',1)

grid on
FS = 'FontSize'; FW = 'FontWeight';
st = sprintf('Eksempel 2 med \\Deltat = 1.58e-5');
title(st,FS,14,FW,'Bold');
xlabel('t',FS,14,FW,'Bold');
ylabel('y_{p}',FS,14,FW,'Bold','Rotation',0);
shg