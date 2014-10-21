%======================= Heunv2 ==========================
% The program computes the velocity of a golf ball
% falling vertically iair and compares with the analytical
% solution. The equation of motion is integrated
% using Heun's method
%=========================================================

clear
nu = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.22  ;  % Density of fluid [kg/m^3]
rol = 1275.0;  % Density of ball [kg/m^3]
d   =  0.041 ; % Diameter of ball [m]
tend = 10.0  ; % Max. time [s]
g = 9.81    ;  % Gravity [N/kg]

% Using Cd = 0.4;
alfa = 7.0e-3;
k1 = sqrt(g/alfa); k2 = sqrt(alfa*g);
hold on
for dt = [0.5 0.2 0.1]
    nsteps = round(tend/dt) + 1;
    va = zeros(nsteps,1); t = va; 
    v = va; feil = v;        % allocate space
    for k = 1:nsteps - 1
        t(k+1) = k*dt;
        vp = v(k) + dt*(g - alfa*v(k)^2);
        v(k+1) = v(k) + 0.5*dt*(2*g - alfa*(v(k)^2 + vp^2));
    end
    % Analytical solution
    va = k1*tanh(k2*t);
   % Compute the relative error
   for k = 2:nsteps
       relerr(k) = abs((va(k) -v(k))/va(k));
   end
   relerr = relerr*100;
   plot(t,relerr)
end
shg
hold off
