%===================== rkule ============================
% The program computes the velocity and the displacement
% of a sphere falling vertically in a fluid. The equation 
% of motion is integrated using RK4C
%========================================================
clear;
clear global A B C d nu Re;
global A B C d nu Re;
n = 2; y = zeros(n,1);
nu = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.22  ;  % Density of fluid [kg/m^3]
rol = 7850.0;  % Density of sphere [kg/m^3]
d   =  0.01 ;  % Diameter of sphere [m]
dt  = 0.1   ;  % Timestep [s]
tend = 10.0  ;  % Max. time [s]
fprintf('       Kinematical viscosity . nu   = %10.3e m^2/s \n',nu );
fprintf('       Density of fluid ...... rof  = %10.3e kg/m^3 \n',rof);
fprintf('       Density of sphere ..... rol  = %10.3e kg/m^3 \n',rol);
fprintf('       Diameter of sphere .... d    = %10.3e m \n',d);
fprintf('       Timestep .............. dt   = %10.3e s \n',dt);
fprintf('       Max. time ............. tend = %10.3e s \n\n',tend);

g = 9.81    ;  % Gravity [N/kg]
ro = rof/rol;
A = 1.0 + 0.5*ro ;
B = (1.0 - ro)*g ;
C = 0.75*ro/d;
nsteps = round(tend/dt);
v = 0.0 ; t = 0.0;Re = 0.0;
fprintf('        t(s)       z(m)         v(m/s)          Re \n\n');
y(1) = 0.0; y(2) = 0.0;
for k = 1:nsteps
   t = k*dt;
   y = rk4c(@fcnkule,t,y,dt); 
   fprintf(' %10.2f  %13.3e %13.3e  %13.3e \n',t,y(1),y(2),Re);
end

   


