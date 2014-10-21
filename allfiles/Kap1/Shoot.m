%===================== SHOOT ===================
% The program computes the movement of a spherical
% projectile in two dimensions. The equations 
% of motion are integrated using RK4C
%======================================================
clear;
global g C d nu Re vfx vfy
n = 4; dy = zeros(n,1); yv = dy;
nu  = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.22   ;  % Density of fluid [kg/m^3]
ros = 7850.0 ;  % Density of sphere [kg/m^3]
d   = 0.01   ;  % Diameter of sphere [m]
dt  = 0.1    ;  % Timestep [s]
v0  = 50.0   ;  % Initial velocity [m/s]
vfx = 0.0    ;  % x-comp. of fluid velocity
vfy = 0.0    ;  % y-comp. of fluid velocity
alf = 45.0   ;  % Angle of elevation (deg.)

fprintf('        Kinematical viscosity . nu   = %10.3e m^2/s \n',nu );
fprintf('        Density of fluid ...... rof  = %10.3e kg/m^3 \n',rof);
fprintf('        Density of sphere ..... ros  = %10.3e kg/m^3 \n',ros);
fprintf('        Diameter of sphere .... d    = %10.3e m \n',d);
fprintf('        Timestep .............. dt   = %10.3e s \n',dt);
fprintf('        Initial velocity....... v0   = %10.3e m/s \n',v0);
fprintf('        x-comp. of fluid vel... vfx  = %10.3e m/s \n',vfx);
fprintf('        y-comp. of fluid vel... vfy  = %10.3e m/s \n',vfy);
fprintf('        Angle of elevation..... alf  = %10.3e deg. \n',alf);

g = 9.81    ;  % Gravity [N/kg]
radf = pi/180;
ro = rof/ros;
C = 0.75*ro/d;
t = 0.0;
vx = v0*cos(alf*radf);
vy = v0*sin(alf*radf);
yv(1) = 0.0; yv(2) = 0.0; yv(3) = vx; yv(4) = vy;
dy = kule2d(t,yv);
fprintf('\n   t(s)      x(m)       y(m)      vx(m/s)      vy(m/s)     Re \n\n');
fs = ' %6.3f %11.3e %11.3e %11.3e %11.3e %11.3e \n';
fprintf(fs,t,yv(1),yv(2),yv(3),yv(4),Re);
%
% === Solving ===
%
k = 0; test = 1;
while ( test > 0.0) 
   k = k + 1;
   t = k*dt;  
   yv = RK4C('kule2d',t,yv,dt);
   test = yv(2);
   if ( test > 0.0 )
      fprintf(fs,t,yv(1),yv(2),yv(3),yv(4),Re);  
   end
end
%
% Interpolate to find the t-value when
% the projectile hits the ground.
%
y2 = yv(2); vy2 = yv(4);
yv = RK4C('kule2d',t,yv,-dt); % One timestep back
dtau = yv(2)*dt/(yv(2) - y2); % Linear estimate
a = (vy2 - yv(4))*0.5/dt;
dtau = -(yv(2) + a*dtau^2)/yv(4); % Improved estimate
yv = RK4C('kule2d',t,yv,dtau);
t = t - dt + dtau ;
fprintf(fs,t,yv(1),yv(2),yv(3),yv(4),Re); % Final step

   


