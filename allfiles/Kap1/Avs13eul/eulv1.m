%===================== Eulv1 ==========================
% The program computes the velocity of a golf ball
% falling vertically in air. The equation set(0,'DefaultTextFontSize',FS);
% of motion is integrated using the basic Euler- method
%=======================================================
clear all; close all; clc;
FS = 20;set(0,'DefaultLineLineWidth',3,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);

nu = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.22  ;  % Density of air [kg/m^3]
rol = 1275.0;  % Density of sphere [kg/m^3]
d   =  0.041 ; % Diameter of ball [m]
tend = 2.0  ; % Max. time [s]
g = 9.81    ;  % Gravity [N/kg]

% Using Cd = 0.4;
alfa = 7.0e-3;
tend = 1.9;
TaylorMax=2
dt = 0.5; % Timestep [s]

nsteps = round(tend/dt) + 1;
v = zeros(nsteps,1); t = v; % allocate space
v(1)= 0.0 ; t(1) = 0.0;
fprintf('       t(s)       v(m/s)        Re \n\n');

for k = 1:nsteps - 1
   t(k+1) = k*dt;
   v(k+1) = v(k) + dt*(g - alfa*v(k)^2);
end

fprintf(' %10.2f  %10.3f %15.3e \n',t,v);

% Analytical solution 
k1 = sqrt(g/alfa); k2 = sqrt(alfa*g);
va = k1*tanh(k2*t);

% Taylor solution
vt=g.*t.*(1-alfa*t.^2+2*alfa^2*g^2*t.^4/15);

if (tend<TaylorMax)
    plot(t,v,t,va,t,vt);
    hl=legend('v','v{_a}','v_t','Location','SouthEast'); 
else
    plot(t,v,t,va);
    hl=legend('v','v{_a}','Location','SouthEast');
end
set(hl,'box','off'); 
grid on

% Title and labels
title(sprintf('Eulers metode med \\Deltat = %4.2f',dt));
xlabel('t(s)'); ylabel('velocity (m/s)');
