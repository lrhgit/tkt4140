% Program Pendel
% Løser pendeligningen på dimensjonsløs form :
% theta'' + sin(theta) = 0
% med startbetingelser :
% theta = theta0 for t = 0
% theta' = 0 for t = 0
% Bruker ode45
clear
% Regner for theta0 = 60
% Svingetiden T for en kvart periode er gitt
% T = K(k) der K er det fullstendige elliptiske
% integtalet av første slag. k = sin(theta0/2)
% for theta0 = 60 blir k = 1/2.
k = 0.5;
T = ellipke(k^2);
y0 = [pi/3; 0]; % Startverdier
tintervall = [0 T];
[t,y] = ode45(@fcn,tintervall,y0);

% Plotter z og v som funksjon av t
plot(t,y(:,1),t,y(:,2),'-.');
grid on

% Tabeller for z og v opptil t = 1.6
options = odeset('RelTol',1.0e-5);
y0 = [pi/3; 0]; % Startverdier
tintervall = [0 :0.1 : 1.6];
[t,y] = ode45(@fcn,tintervall,y0,options);
fprintf(' %5.1f  %13.4e  %13.4e \n',[t y]');

