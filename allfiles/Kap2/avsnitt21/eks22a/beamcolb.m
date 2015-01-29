%============== Program beamcolv2 =========================
% The program computes the deflection of a beam-column
% with a uniform lateral load of intensity q and an axial
% compressive force P using  a shooting technique.
% The differential equation is solved using ode45.
% In this version, the length of the beam is 2 with
% -1 <= x <= 1.
% Boundary conditions : u(-1) = u(1) = 0
%                 or  : u'(0) = 0 , u(1) = 0 
% s = u(0)
% Guessing s0 = 0.5, s1 = 1.2
%
% Using function fcnbeam2
%
%========================================================
clear; clear global theta2;
global theta2 
% === Initialize ===
theta = 1.0   ;  % Load-parameter
% Note : Buckling for theta = Pi/2
theta2 = theta^2;
fprintf('  Load-parameter theta  = %10.3e \n',theta);
x0 = 0; x1 = 1;
xspan = [x0 x1];
s0 = 0.2;
s1 = 0.5;
options = odeset('RelTol',1.0e-5);
% === Compute fi0;
u0 = [s0 0];
[x,u] = ode45(@fcnbeam2,xspan, u0,options);
fi0 = u(end,1);
% === Compute fi1;
u0 = [s1 0];
[x,u] = ode45(@fcnbeam2,xspan, u0,options);
fi1 = u(end,1);
s = (fi1*s0 - fi0*s1)/(fi1 - fi0);
% === Compute a table for u and du/dx;
xspan = (x0: 0.1:x1);
u0 = [s 0];
[x,u] = ode45(@fcnbeam2,xspan, u0,options);
%x = (0 : dx : 1.0)';
% --- Analytical solution ua 
ua = (cos(theta*x)/cos(theta) -1)/theta2 - (1-x.^2)/2;

% === Output ===
fprintf('\n    x       u-comput.      u-analyt.   \n\n'); 
fprintf(' %7.3f   %12.4e  %12.4e \n',[x u(:,1) ua]'); 

% === Plotting
xspan = [x0 x1];
u0 = [s 0];
[x,u] = ode45(@fcnbeam2,xspan, u0,options);
plot(x,u(:,1),x,u(:,2),'-.');
grid



