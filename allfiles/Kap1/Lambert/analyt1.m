% Program analyt1
% Plotter den analytiske løsningen av 
% oppgave 1.6.1 i Lambert.
% 
clear; close;
x = linspace(0,0.2);
v1 = exp(-40*x);
v2 = cos(40*x);
v3 = sin(40*x);
y1 = (exp(-2*x) + v1.*(v2 - v3))/2;
y2 = (exp(-2*x) + v1.*(-v2 + v3))/2;
y3 = -v1.*(v2 + v3);
plot(x,y1,x,y2,x,y3);
grid
y1 = y1'; y2 = y2'; y3 = y3';
[y1 y2 y3]
% Tables of z and v 
% options = odeset('RelTol',1.0e-5);
% y0 = [1; 0]; % Initial values
% tspan = (0 : 1 : 10);
% [t,y] = ode45(@fcn2,tspan,y0,options);
% fprintf(' %5.1f  %13.4e  %13.4e \n',[t y]');

