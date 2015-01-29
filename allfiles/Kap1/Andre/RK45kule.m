%===================== RK45KULE ===================
% The program comutes the velocity of a sphere
% falling vertically in a fluid. The equation 
% of motion is integrated using ODE45
%======================================================
clear;
global A B C d nu Re
nu = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.22  ;  % Density of fluid [kg/m^3]
rol = 7850.0;  % Density of sphere [kg/m^3]
d   =  0.01 ;  % Diameter of sphere [m]
dt  = 0.05   ;  % Timestep [s]
tend = 2.0  ;  % Max. time [s]
options = odeset('RelTol',1.0e-5);
tspan = [0.0 : 0.1 :tend];
y0 =[ 0.0 0.0];
fprintf(' Kinematical viscosity . nu   = %10.3e m^2/s \n',nu );
fprintf(' Density of fluid ...... rof  = %10.3e kg/m^3 \n',rof);
fprintf(' Density of sphere ..... rol  = %10.3e kg/m^3 \n',rol);
fprintf(' Diameter of sphere .... d    = %10.3e m \n',d);
fprintf(' Timestep .............. dt   = %10.3e s \n',dt);
fprintf(' Max. time ............. tend = %10.3e s \n\n',tend);

g = 9.81    ;  % Gravity [N/kg]
ro = rof/rol;
A = 1.0 + 0.5*ro ;
B = (1.0 - ro)*g ;
C = 0.75*ro/d;
[t,y] = ode45('fcn',tspan,y0,options);
fprintf('       t(s)       v(m/s)        Re \n\n');
fprintf(' %10.2f  %13.3e %13.3e \n',[t,y(:,1),y(:,2)]');
end

   


