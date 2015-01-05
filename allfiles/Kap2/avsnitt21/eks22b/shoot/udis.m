%================== Program udis =========================
% The program computes the deflection of a beam-column
% with a uniform lateral load of intensity q and an axial
% compressive force P using  a shooting technique.
% The bending stiffiness EI varies according to
% I = I0/(1+x^2), -1 <= x <= 1 , E is constant
% The differential equation is solved using ode45.
% In this version, the length of the beam is 2 with
% -1 <= x <= 1.
% Boundary conditions :u(-1) = u(1) = 0
%                 or  : u'(0) = 0 , u(1) = 0 
% s = m(0)
% Guessing s0 = 0.5, s1 = 1.2
%
% Using function fcnbeamu
%
%========================================================
clear; 
%
% === Initialize ===
x0 = 0; x1 = 1;
xspan = [x0 x1];
s0 = 0.5;
s1 = 1.2;
options = odeset('RelTol',1.0e-5);
% === Compute fi0;
u0 = [s0 0];
[x,u] = ode45(@fcnbeamu,xspan, u0,options);
fi0 = u(end,1);
% === Compute fi1;
u0 = [s1 0];
[x,u] = ode45(@fcnbeamu,xspan, u0,options);
fi1 = u(end,1);
s = (fi1*s0 - fi0*s1)/(fi1 - fi0);
% === Compute a table for u  ;
xspan = (x0: 0.1:x1);
u0 = [s 0];
[x,u] = ode45(@fcnbeamu,xspan, u0,options);

%=== Output ===
fprintf('\n    x       u-comput.     \n\n'); 
fprintf(' %7.3f   %12.4e   \n',[x u(:,1) ]'); 

% === Plotting
% xspan = [x0 x1];
% m0 = [s 0];
% [x,m] = ode45(@fcnbeam,xspan, m0,options);
% plot(x,m(:,1),x,u(:,2),'-.');
% grid



