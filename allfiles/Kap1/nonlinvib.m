function nonlinvib
% Large amplitude vibration of a mass on a frictionless
% foundation together with the analytical solution for the
% small-amplitude case.
% Equation: u''(t) + 2*u(t)*(fac -1)/fac = 0
% with fac = sqrt(1 + u(t)^2)
% Initial conditions: u(0) = -u0, u'(0) = 0
% Small-amplitude equation:
%     u''(t) + u(t)^3 = 0
% 
% Integration using ODE45
%==========================================================
clear;
y0 = zeros(2,1); y = y0; % Initialize
u0 = input('u0  = '); 
ystart = -u0; % Initial displacement
tend = 40.0;  % Max. time  
tspan = [0 tend];
y0 =[ -u0 0.0]';
fprintf(' Initial displacement ...... ystart = %10.3e \n',ystart);
fprintf(' Max. time ................. tend   = %10.3e \n',tend);
options = odeset('RelTol',1.0e-5,'Events',@event1);
[t,y,te,ye,ie] = ode45(@fcn,tspan,y0,options);
fprintf('\n The mass is first time in origo when t = %6.4f',te);
fprintf('\n The velocity in origo is %6.4f',ye(2));
p = 4*te;
fprintf('\n The period = %6.4f \n',p);
% === Small-amplitude solution
ta = u0*t;
[sn,cn,dn] = ellipj(ta,0.5); % Using Jacobi elliptical functions
ua = -u0*cn; % Analytical solution
% === Plotting ===
FS = 'FontSize'; FW = 'FontWeight';
plot(t,y(:,1),'k',t,ua,'k-.');
grid 
xlabel('t',FS,14,FW,'Bold')
ylabel('u , u{_a}',FS,14,FW,'Bold')
st = sprintf('Nonlinear Vibration.  Period = %4.3f',p);
title(st,FS,14)
legend('u','u{_a}',0)
% ----------------------------------------------------
function dydt = fcn(t,y)
dydt = zeros(size(y));
fac = sqrt(1 + y(1)^2);
dydt(1) = y(2);
dydt(2) = - 2*y(1)*(fac - 1)/fac;
%-------------------------------------------------------
function [value, isterminal,direction] = event1(t,y)
value = y(1);
isterminal = 1;
direction = 0;

   


