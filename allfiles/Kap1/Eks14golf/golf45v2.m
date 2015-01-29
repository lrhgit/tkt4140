function golf45v2
%=================================================================
% The program computes the trajectory of a golf ball.
% The equations of motion are integrated using ODE45
% The program uses the functions fcngball45,cdcldata and events
% Dragdata from function cdcldata.
% ========================= Reference ====================== 
% P. W. Bearman & J. K. Harvey :
% "Golf Ball Aerodynamics",
% Aeronautical Quarterly vol. 27, 1976, pp. 112-122
% ==========================================================
clear; clear global g C vfx vfy nrpm;
global g C vfx vfy nrpm ;
neq = 4; y = zeros(neq,1); 
nu  = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.20   ;  % Density of fluid [kg/m^3]
ros = 1260.0 ;  % Density of sphere [kg/m^3]. m = 46g
d   = 0.041   ;  % Diameter of sphere [m]
v0  = 61.0   ;  % Initial velocity [m/s]
vfx = 0.0    ;  % x-comp. of fluid velocity
vfy = 0.0    ;  % y-comp. of fluid velocity
nrpm = 6000  ;  % Number of revolutions pr. minute

fprintf('        Kinematical viscosity . nu   = %10.3e m^2/s \n',nu );
fprintf('        Density of fluid ...... rof  = %10.3e kg/m^3 \n',rof);
fprintf('        Density of sphere ..... ros  = %10.3e kg/m^3 \n',ros);
fprintf('        Diameter of sphere .... d    = %10.3e m \n',d);
fprintf('        Initial velocity....... v0   = %10.3e m/s \n',v0);
fprintf('        x-comp. of fluid vel... vfx  = %10.3e m/s \n',vfx);
fprintf('        y-comp. of fluid vel... vfy  = %10.3e m/s \n',vfy);
fprintf('        Number of revolutions. nrpm  = %10.3e rpm \n',nrpm);

g = 9.81    ;  % Gravity [N/kg]
radf = pi/180;
ro = rof/ros;
C = 0.75*ro/d;
alf = 15; % Angle of elevation
vx = v0*cos(alf*radf);
vy = v0*sin(alf*radf);
tint = [0 10]; % timeinterval
y0 = [0.0 ; 0.0; vx; vy]; % Initial values
options = odeset('RelTol',1.0e-5,'Refine',8,'Events',@events);
[t,y,te,ye,ie] = ode45(@fcngball45,tint,y0,options);
% === Plotting the trajectories ===
FS = 'FontSize'; FW = 'FontWeight';
st = sprintf('Utgangshastighet = %5.1f m/s . Spinn = %5.0f o/min',v0,nrpm);
plot(y(:,1),y(:,2),'k');
daspect([1 1 1]);
axis([0 220 0 70])
xlabel('x(m)',FS,12,FW,'Bold')
ylabel('y(m)',FS,12,FW,'Bold')
title(st,FS,11,FW,'Bold')
hold on
for alf = [20 25 30] % Angle of elevation
    vx = v0*cos(alf*radf);
    vy = v0*sin(alf*radf);
    y0 = [0.0 ; 0.0; vx; vy]; % Initial values
    [t,y,te,ye,ie] = ode45(@fcngball45,tint,y0,options);
    plot(y(:,1),y(:,2),'k');
end
hold off
%========== fcngball45 =========
function dydt = fcngball45(t,y)
% Called by ode45
% Dragdata from function cdcldata
dydt = zeros(size(y));
global g C vfx vfy nrpm;
vrx = y(3) - vfx;
vry = y(4) - vfy;
vr = sqrt(vrx^2 + vry^2);
[cd,cl] = cdcldata(vr,nrpm);
f1 = C*vr*cd; f2 = C*vr*cl;
dydt(1) = y(3);
dydt(2) = y(4);
dydt(3) = - f1*vrx - f2*vry;
dydt(4) = - g + f2*vrx - f1*vry;
%================ events ===========================
function [value,isterminal,direction] = events(t,y)
value = y;
isterminal = [0; 1; 0 ;0];
direction =  [0; -1; 0; 0];

