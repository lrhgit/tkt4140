%===================== EKULE ============================
% The program computes the velocity and displacement
% of a sphere falling vertically in a fluid. The equations 
% of motion are integrated using the forward Euler-method
% The first part output a table from t = 0 to t = tend.
% The second part integrates until the terminal velocity
% is attained and plots v and z versus t.
% We have retained the absolute sign for v even if
% it not needed in this case.
% Drag computed from function CDkule
%=========================================================
clear all; close all;
clc;

nu = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.22  ;  % Density of fluid [kg/m^3]
rol = 1275.0;  % Density of sphere [kg/m^3]
d   = 0.041 ;  % Diameter of sphere [m] (golf ball)
dt  = 0.1   ;  % Timestep [s]
tend = 2.0  ;  % Max. time [s] for table
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
nsteps = round(tend/dt);
z = 0.0 ;v = 0.0 ; t = 0.0;
fprintf('       t(s)       v(m/s)      z(m)        Re \n\n');

%%  ===== OUTPUT A TABLE ====
for k = 1:nsteps
   t = k*dt;
   va = abs(v); 
   Re = va*d/nu;
   CD = CDkule(Re);
   f = (B - C*v*va*CD)/A;
   z = z + dt*v;
   v = v + dt*f;
   fprintf(' %10.2f  %10.3f %10.3f %15.3e \n',t,v,z,Re);
end

%% === COMPUTE UNTIL THE TERMINAL VELOCITY IS ATTAINED ====
% We collect v, z and t in vectors for plotting
v(1) = 0.0; z(1) = 0.0; t(1) = 0.0; vt = 1.0;
epsi = 5.0e-3; stest = 1; k = 0;
while stest > epsi
    k = k + 1;
    t(k+1) = k*dt;
    va = abs(v(k)); 
    Re = va*d/nu;
    CD = CDkule(Re);
    f = (B - C*v(k)*va*CD)/A;
    z(k+1) = z(k) + dt*v(k);
    v(k+1) = v(k) + dt*f;
    if k > 1
        vt = sqrt(B/(CD*C));
    else
        vt = 0;
    end
    stest = abs((vt - v(k+1))/v(k+1));
end

%% Report results
fprintf('\n Terminal velocity = %7.3f (m/s) at t = %7.3f s  \n',v(end),t(end));
%plot(t,v,'k',t,z,'k');
h=plot(t,v);
set(h,'linewidth',2);
FS = 'FontSize'; FW = 'FontWeight';
xlabel('t(s)',FS,14)
ylabel('v(m/s)',FS,14) 
st= sprintf('Falling sphere. Euler''s method with \\Deltat = %4.2f',dt);
title(st,FS,14)


