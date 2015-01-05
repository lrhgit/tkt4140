function sball45v2
%=================================================================
% The program computes the trajectory of a smooth ball with no lift
% The equations of motion are integrated using ODE45
% The program uses the functions fcnsball45,CDkule and events.
% This version computes and plots trajectories for several
% values of the elevation-angle. Dragdata from function CDkule
%==================================================================
clear all; close all;
clc;

global g C d nu vfx vfy
neq = 4; y = zeros(neq,1); 
nu  = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.20   ;  % Density of fluid [kg/m^3]
ros = 1275.0 ;  % Density of sphere [kg/m^3]. Mass m = 50g
d   = 0.04   ;  % Diameter of sphere [m]
v0  = 50.0   ;  % Initial velocity [m/s]
vfx = 0.0    ;  % x-comp. of fluid velocity
vfy = 0.0    ;  % y-comp. of fluid velocity

fprintf('        Kinematical viscosity . nu   = %10.3e m^2/s \n',nu );
fprintf('        Density of fluid ...... rof  = %10.3e kg/m^3 \n',rof);
fprintf('        Density of sphere ..... ros  = %10.3e kg/m^3 \n',ros);
fprintf('        Diameter of sphere .... d    = %10.3e m \n',d);
fprintf('        Initial velocity....... v0   = %10.3e m/s \n',v0);
fprintf('        x-comp. of fluid vel... vfx  = %10.3e m/s \n',vfx);
fprintf('        y-comp. of fluid vel... vfy  = %10.3e m/s \n',vfy);

g = 9.81    ;  % Gravity [N/kg]
radf = pi/180;
ro = rof/ros;
C = 0.75*ro/d;
alfa = 20; % Angle of elevation
vx = v0*cos(alfa*radf);
vy = v0*sin(alfa*radf);
tint = [0 10]; % timeinterval
y0 = [0.0 ; 0.0; vx; vy]; % Initial values
options = odeset('RelTol',1.0e-5,'Refine',8,'Events',@events);
[t,y,te,ye,ie] = ode45(@fcnsball45,tint,y0,options);


%% === Plotting the trajectories ===
FS = 14; FW = 'FontWeight';
h(1,:)=plot(y(:,1),y(:,2));
nplt=1; %plot counter

% daspect([1 1 1]);
xlim([0 120])
ylim([0 40])
% set(gca,'FontName','Century Schoolbook');
xlabel('x(m)','FontSize',FS,FW,'Bold')
ylabel('y(m)','FontSize',FS,FW,'Bold')
title('Utgangshastighet = 50m/s','FontSize',FS,FW,'Bold')
hold on

for alfa = [30 40 45] % Angle of elevation
    vx = v0*cos(alfa*radf);
    vy = v0*sin(alfa*radf);
    y0 = [0.0 ; 0.0; vx; vy]; % Initial values
    [t,y,te,ye,ie] = ode45(@fcnsball45,tint,y0,options);
    nplt=nplt+1;
    h(nplt,:)=plot(y(:,1),y(:,2));
end
hold off
set(h(:,:),'linewidth',2);


function dydt = fcnsball45(t,y)
% Used by the program sball45v2
% Called by ode45
% Dragdata from function CDkule
dydt = zeros(size(y));
global g C d nu vfx vfy;
vrx = y(3) - vfx;
vry = y(4) - vfy;
vr = sqrt(vrx^2 + vry^2);
Re = vr*d/nu;
CD = CDkule(Re);
f = C*vr*CD;
dydt(1) = y(3);
dydt(2) = y(4);
dydt(3) = - f*vrx;
dydt(4) = - g -f*vry;

function [value,isterminal,direction] = events(t,y)
value = y;
isterminal = [0; 1; 0 ;0];
direction =  [0; -1; 0; 0];

