% Program versjon1 
% Bruker ode45
clear
y0 = [8; 0]; % Startverdier
tintervall = [0 70.0];
[t,y] = ode45('fcn',tintervall,y0);

% Plotter z og v = dz/dt som funksjon av t
plot(t,y(:,1),t,y(:,2),'-.')
