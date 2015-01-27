%===================== Heunz2 ==========================
% The program computes the velocity and displacement of 
% a golf ball falling vertically in air and compares 
% with the analytical solution. The equation of motion 
% is integrated using Heun's method.
%=======================================================

clear
nu = 1.5e-5 ;   % Kinematical viscosity of air [m^2/s]
rof = 1.22  ;   % Density of air [kg/m^3]
rol = 1275.0;   % Density of golf ball [kg/m^3]
d   =  0.041 ;  % Diameter of ball [m]
tend = 10.0  ;  % Max. time [s]
g = 9.81    ;   % Gravity [N/kg]

% Using Cd =  0.4
alfa = 7.0e-3;
k1 = sqrt(g/alfa); k2 = sqrt(alfa*g);
hold on
for dt = [ 0.1 0.2 0.5 ]
    nsteps = round(tend/dt) + 1;
    % --- Allocate space for vectors
    za = zeros(nsteps,1); t = za; z = za; 
    v = za; relerr = v;
    % --- Heun's method
    for k = 1:nsteps - 1
        t(k+1) = k*dt;
        % --- Predictor
        vp = v(k) + dt*(g - alfa*v(k)^2);
        % --- Corrector
        z(k+1) = z(k) + 0.5*dt*(v(k) + vp);
        v(k+1) = v(k) + 0.5*dt*(2*g - alfa*(v(k)^2 + vp^2));
    end
    % Analytical solution
    za = log(cosh(k2*t))/alfa;
   % Compute relative error 
   for k = 2:nsteps
       relerr(k) = abs((za(k) - z(k))/za(k));
       %relerr(k) = abs(za(k) - z(k)); % abs. err
   end
   relerr = relerr*100; % Error in percent
   plot(t,relerr)
end
shg
hold off
