% Program lambert1
% Løser oppgave 1.6.1 i Lambert numerisk
% Using ode45
clear; close;
y0 = [1; 0; -1]; % Initial values
xspan = [0 1.0];
options = odeset('RelTol',1.0e-5);
[x,y] = ode45(@fcnlam1,xspan,y0);
plot(x,y(:,1),x,y(:,2),x,y(:,3));
grid

% Tables of z and v 
% options = odeset('RelTol',1.0e-5);
% y0 = [1; 0]; % Initial values
% tspan = (0 : 1 : 10);
% [t,y] = ode45(@fcn2,tspan,y0,options);
% fprintf(' %5.1f  %13.4e  %13.4e \n',[t y]');

