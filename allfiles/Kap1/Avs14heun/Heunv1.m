%===================== Heunv1 ==========================
% The program computes the velocity of a golf ball
% falling vertically in air. The equation 
% of motion is integrated using Heun's method
%=======================================================

clear
nu = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.22  ;  % Density of air [kg/m^3]
rol = 1275.0;  % Density of ball [kg/m^3]
d   =  0.041 ; % Diameter of ball [m]
tend = 10.0  ; % Max. time [s]
g = 9.81    ;  % Gravity [N/kg]

% Using Cd = 0.4;
alfa = 7.0e-3;
tend = 10;
dt = 0.5;
nsteps = round(tend/dt) + 1;
v = zeros(nsteps,1); t = v; feil = v; % allocate space
v(1)= 0.0 ; t(1) = 0.0; feil(1) = 0;
%fprintf('       t(s)       v(m/s)        Re \n\n');
for k = 1:nsteps - 1
   t(k+1) = k*dt;
   vp = v(k) + dt*(g - alfa*v(k)^2);
   v(k+1) = v(k) + 0.5*dt*(2*g - alfa*(v(k)^2 + vp^2));
end
%fprintf(' %10.2f  %10.3f %15.3e \n',t,v,Re);
% Analytical solution
k1 = sqrt(g/alfa); k2 = sqrt(alfa*g);
va = k1*tanh(k2*t);
plot(t,v,t,va);
shg
% Compute the relative error
for k = 2:nsteps
   relerr(k) = abs((va(k) -v(k))/va(k));
end
relerr = relerr*100;
plot(t,relerr)
