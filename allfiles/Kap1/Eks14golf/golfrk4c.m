%========================= GOLFRK4C ========================
% The program computes the trajectory of a golf ball.
% The equations of motion are integrated using RK4C
% The system of equations is programmed in function fcngball
% Dragdata from function cdcldata.
% Output: Tables of coordinates ,velocities
% and Reynolds number as function of time. (No plotting)
% ========================= Reference ====================== 
% P. W. Bearman & J. K. Harvey :
% "Golf Ball Aerodynamics",
% Aeronautical Quarterly vol. 27, 1976, pp. 112-122
% ==========================================================
clear;
global g C d nu Re vfx vfy nrpm
n = 4; yvec = zeros(n,1); ystart = yvec;
nu  = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.20   ;  % Density of fluid [kg/m^3]
ros = 1275.0 ;  % Density of ball [kg/m^3] . m = 50g
d   = 0.041   ;  % Diameter of ball [m]
dt  = 0.1    ;  % Timestep [s]
v0  = 61.0   ;  % Initial velocity [m/s]
vfx = 0.0    ;  % x-comp. of fluid velocity
vfy = 0.0    ;  % y-comp. of fluid velocity
alf = 21.0    ;  % Angle of elevation (deg.)
nrpm = 3500.0 ;  % number of revolutions pr. minute

fprintf('        Kinematical viscosity . nu   = %10.3e m^2/s \n',nu );
fprintf('        Density of fluid ...... rof  = %10.3e kg/m^3 \n',rof);
fprintf('        Density of ball ....... ros  = %10.3e kg/m^3 \n',ros);
fprintf('        Diameter of ball ...... d    = %10.3e m \n',d);
fprintf('        Timestep .............. dt   = %10.3e s \n',dt);
fprintf('        Initial velocity....... v0   = %10.3e m/s \n',v0);
fprintf('        x-comp. of fluid vel... vfx  = %10.3e m/s \n',vfx);
fprintf('        y-comp. of fluid vel... vfy  = %10.3e m/s \n',vfy);
fprintf('        Angle of elevation..... alf  = %10.3e deg. \n',alf);
fprintf('        Number of revolutions.. nrpm = %10.3e rpm. \n',nrpm);

g = 9.81    ;  % Gravity [N/kg]
radf = pi/180;
ro = rof/ros;
C = 0.75*ro/d;
t = 0.0;
vx = v0*cos(alf*radf);
vy = v0*sin(alf*radf);
ystart(1) = 0.0; ystart(2) = 0.0; ystart(3) = vx; ystart(4) = vy;
fprintf('\n   t(s)      x(m)       y(m)      vx(m/s)      vy(m/s)     Re \n\n');
fs = '%6.3f %11.3e %11.3e %11.3e %11.3e %11.3e \n';
fprintf(fs,t,ystart(1),ystart(2),ystart(3),ystart(4),Re);
%
% === Solving ===
%
k = 0; test = 1;
while ( test > 0.0) 
   k = k + 1;
   t = k*dt;  
   yvec = RK4C('fcngball',t,ystart,dt);
   ystart = yvec; % Starting values for the next step
   test = yvec(2); % Testing for positive y-value
   if ( test > 0.0 )
      fprintf(fs,t,yvec(1),yvec(2),yvec(3),yvec(4),Re);  
   end
end
%
% Interpolate to find the t-value when
% the ball hits the ground.
%
vy2 = yvec(4);
yvec = RK4C('fcngball',t,yvec,-dt); % Back one timestep
tau0l = -yvec(2)/yvec(4); % Linear estimate
a = (vy2 - yvec(4))/dt;
tau0 = tau0l*(1 - 0.5*a*tau0l/yvec(4)); % Improved estimate
yvec = RK4C('fcngball',t,yvec,tau0); % Final step
t = t - dt + tau0 ;
fprintf(fs,t,yvec(1),yvec(2),yvec(3),yvec(4),Re); 

   


