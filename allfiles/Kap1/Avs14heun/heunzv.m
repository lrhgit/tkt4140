%================== Program heunzv ======================
% The program computes the velocity and the displacement
% of a golf ball falling vertically in air. The equation 
% of motion is integrated using Heun's method
%========================================================
clear all; close all; clc;
FS = 20; set(0,'DefaultLineLineWidth',3,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);
% Using Cd = 0.4;
g = 9.81; % Gravity [N/kg]
tend = 10;
alfa = 7.0e-3;
dt = 1.5; % Timestep
steps = round(tend/dt) + 1;
v = zeros(steps,1); t = v; z = v; % allocate space
v(1)= 0.0 ; t(1) = 0.0; z(1) = 0;
for n = 1:steps - 1
   t(n+1) = n*dt;
   % --- Predictor
   % zp = z(n) + dt*v(n) not used
   vp = v(n) + dt*(g - alfa*v(n)^2);
   % --- Corrector
   z(n+1) = z(n) + 0.5*dt*(v(n) + vp);
   v(n+1) = v(n) + 0.5*dt*(2*g - alfa*(v(n)^2 + vp^2));
end
% Analytical solution
k1 = sqrt(g/alfa); k2 = sqrt(alfa*g);
va = k1*tanh(k2*t);
za = log(cosh(k2*t))/alfa;
fprintf('    t(s)    v(m/s)   va(m/s)   z(m)    za(m) \n\n');
fprintf(' %6.2f  %8.3f %8.3f %8.2f %8.2f\n',[t v va z za]');
plot(t,v,t,va)

