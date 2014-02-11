%================== Program momentdis =========================
% The program computes the moment-distribution of a beam-column
% with a uniform lateral load of intensity q and an axial
% compressive force P using  a shooting technique.
% The bending stiffiness EI varies according to
% I = I0/(1+x^n), -1 <= x <= 1 , E is constant, n = 2,4,6,..
%
% Equation :
%      m''(x) + (1 + x^n)*m(x) = -1 ,  -1 <= x <= 1.
%
% Boundary conditions used:  m'(0) = 0 , m(1) = 0 
%                
% s = m(0)
% Guessing s0 = 0.0, s1 = 1.0
%
% The differential equation is solved using ode45.
% Using function fcnbeam2
%
%==============================================================
clear; clear gobal n;
global n;
%
% === Initialize ===
n = 2;
x0 = 0; x1 = 1;
xspan = [x0 x1];
s0 = 0.0;
s1 = 1.0;
options = odeset('RelTol',1.0e-5);
% === Compute fi0;
m0 = [s0 0];
[x,m] = ode45(@fcnbeam2,xspan, m0,options);
fi0 = m(end,1);
% === Compute fi1;
m0 = [s1 0];
[x,m] = ode45(@fcnbeam2,xspan, m0,options);
fi1 = m(end,1);
s = (fi1*s0 - fi0*s1)/(fi1 - fi0)
% === Compute a table for m and dm/dx;
xspan = (x0: 0.1:x1);
m0 = [s 0];
[x,m] = ode45(@fcnbeam2,xspan, m0,options);
u = m(:,1) -(1-x.^2)/2;

%=== Output ===
fprintf('\n    x        m-comput.        u    \n\n'); 
fprintf(' %7.3f   %12.4e  %12.4e \n',[x m(:,1) u ]'); 

% === Plotting
% xspan = [x0 x1];
% m0 = [s 0];
% [x,m] = ode45(@fcnbeam2,xspan, m0,options);
% plot(x,m(:,1),x,u(:,2),'-.');
% grid



