%====================== Eulz2 ==========================
% The program computes the velocity and displacement of 
% a golf ball falling vertically in air and compares 
% with the analytical solution. The equation of motion 
% is integrated using Euler's method
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
k1 = sqrt(g/alfa); k2 = sqrt(alfa*g);
hold on
for dt = [ 0.1]
    nsteps = round(tend/dt) + 1;
    za = zeros(nsteps,1); t = za; z = za; 
    v = za; feil = v;        % allocate space
    for k = 1:nsteps - 1
        t(k+1) = k*dt;
        z(k+1) = z(k) + dt*v(k);
        v(k+1) = v(k) + dt*(g - alfa*v(k)^2);
    end
    % Analytical solution
    za = log(cosh(k2*t))/alfa;
   
   % Compute the absolute or the relative error
   for k = 1:nsteps
       relerr(k) = abs(za(k) - z(k));
   end
   %relerr = relerr*100;
   plot(t,relerr)
end
shg
hold off
