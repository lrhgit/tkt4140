function golf45v1
%=================================================================
% The program computes the trajectory of a golf ball.
% The equations of motion are integrated using ODE45
% The program uses the functions fcngball45,cdcldata and events
% Drag- and lift-data from function cdcldata.
% ========================= Reference ====================== 
% P. W. Bearman & J. K. Harvey :
% "Golf Ball Aerodynamics",
% Aeronautical Quarterly vol. 27, 1976, pp. 112-122
% ==========================================================
clear;
global g C d nu Re vfx vfy nrpm
neq = 4; y = zeros(neq,1); 
nu  = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.20   ;  % Density of air [kg/m^3]
ros = 1580.0 ;  % Density of ball [kg/m^3]
d   = 0.041  ;  % Diameter of ball [m]
v0  = 50.0   ;  % Initial velocity [m/s]
vfx = 0.0    ;  % x-comp. of air velocity
vfy = 0.0    ;  % y-comp. of air velocity
alf = 11.0   ;  % Angle of elevation (deg.)
nrpm = 3500  ;  % Number of revolutions pr. minute

fprintf('     Kinematical viscosity . nu    = %10.3e m^2/s \n',nu );
fprintf('     Density of fluid ...... rof   = %10.3e kg/m^3 \n',rof);
fprintf('     Density of sphere ..... ros   = %10.3e kg/m^3 \n',ros);
fprintf('     Diameter of sphere .... d     = %10.3e m \n',d);
fprintf('     Initial velocity....... v0    = %10.3e m/s \n',v0);
fprintf('     x-comp. of fluid vel... vfx   = %10.3e m/s \n',vfx);
fprintf('     y-comp. of fluid vel... vfy   = %10.3e m/s \n',vfy);
fprintf('     Angle of elevation..... alf   = %10.3e deg. \n',alf);
fprintf('     Number of revolutions.. nrpm  = %10.3e rpm. \n',nrpm);

g = 9.81    ;  % Gravity [N/kg]
radf = pi/180;
ro = rof/ros;
C = 0.75*ro/d;
vx = v0*cos(alf*radf);
vy = v0*sin(alf*radf);
tint = [0 10]; % Timeinterval
y0 = [0.0 ; 0.0; vx; vy]; % Initial values
options = odeset('RelTol',1.0e-5,'Refine',8,'Events',@events);
[t,y,te,ye,ie] = ode45(@fcngball45,tint,y0,options);
% === Plotting the trajectory ===
FS = 'FontSize'; FW = 'FontWeight';
plot(y(:,1),y(:,2));
%daspect([1 1 1]);
ylim([0 60])
xlim([0 220])
xlabel('x(m)',FS,14,FW,'Bold')
ylabel('y(m)',FS,14,FW,'Bold')
%========= fcngball45 ===========
function dydt = fcngball45(t,y)
% Called by ode45
% Drag- and lift-data from function cdcldata
dydt = zeros(size(y));
global g C d nu Re vfx vfy nrpm;
vrx = y(3) - vfx;
vry = y(4) - vfy;
vr = sqrt(vrx^2 + vry^2);
Re = vr*d/nu;
[cd,cl] = cdcldata(vr, nrpm);
f1 = C*vr*cd; f2 = C*vr*cl;
dydt(1) = y(3);
dydt(2) = y(4);
dydt(3) = - f1*vrx - f2*vry;
dydt(4) = - g + f2*vrx - f1*vry;
% ================ events ========================
function [value,isterminal,direction] = events(t,y)
value = y;
isterminal = [0; 1; 0 ;0];
direction =  [0; -1; 0; 0];

