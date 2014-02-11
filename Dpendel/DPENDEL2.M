 function dpendel2
% 
% This is a special version of DPENDEL
% Here we focus on the influence of changing
% the accuracy in the equation solver.
%
% To write the equations as a system of 4 1.order equations,
% we put theta1 = y1 and theta2 = y3
% giving d(theta1)/dt = d(y1)/dt = y2
% and d(theta2)/dt = d(y3)/dt = y4
%
% Internal functions : FCN, MASS and PLOTTING 

% Programmer : Jan B. Aarseth,
% Institute of Structural Mechanics, NTNU

global g m1 m2 l1 l2;
g = 9.81; % Gravity force

% ==== Data for masses and length of pendulums ===
m1 = 0.5 ; m2 = 2.0;
l1 = 1.5; l2 = 1.0;

% === Time-data ====
tstopp = 35;
antall = 350;
tstart = 0.0; 
tspan = linspace(tstart,tstopp,antall);

% === Initial values for angles ====
% Standard values : theta1 = 45 og theta2 = 175
theta1g = 45;
theta1 = theta1g*pi/180; % theta1(radians)
theta2g = 150;
theta2 = theta2g*pi/180; % theta2 (radians)

% === Initial values of angular velocities ===
dtheta1 = 0.0; 
dtheta2 = 0.0; 

L = l1 + l2; % total length of pendulum
FW = 'FontWeight'; FS = 'FontSize'; LW = 'LineWidth';

% ===== Plotting of initial configuration ====
clf reset;
x1 = l1*sin(theta1); y1 = l1*cos(theta1);
x2 = x1 + l2*sin(theta2); y2 = y1 + l2*cos(theta2);
xl1 = [0; x1]; yl1 = [0 ; y1];
xl2 = [x1; x2]; yl2 = [y1 ; y2];
plot(xl1,yl1,xl2,yl2,'-ro',LW,3);
axis([-L L -L L],'ij');
grid on
xlabel('x',FW,'Bold',FS,14);
ylabel('y','Rotation',0,FW,'Bold',FS,14);
st = sprintf('l_1 = %3.2f , l_2 = %3.2f ,  m_1 = %3.2f ,  m_2 = %3.2f ',l1,l2,m1,m2);
title(st,FS,13);
st1 = sprintf('\\theta_1 = %3.2f\\circ ',theta1g);
xt = -L*0.9; yt = -L*0.8;
text(xt,yt,st1,FS,14,FW,'Bold');
st2 = sprintf('\\theta_2 = %3.2f\\circ ',theta2g);
yt = -L*0.6;
text(xt,yt,st2,FS,13,FW,'Bold');
st3 = ('Press a key ....');
yt = L*0.9;
text(xt,yt,st3,'Color',[0.5 0 1],FS,13);
pause;

%==== Solving the differential equations ===
clf reset;
epsrel = 1.0e-5;
fac = 180/pi;
h = waitbar(0,'Computing ....');
antall = 4; deler = 4/2;
for k = 1:antall
    epsrel = epsrel*0.1;
    y0 =[theta1; dtheta1; theta2; dtheta2];
    options = odeset('RelTol',epsrel,'Mass',@mass);
    [t,y] = ode113(@fcn,tspan,y0,options);
    % === Plotting the angles as a function of time ===
    subplot(deler,deler,k)
    st = sprintf('RelTol = %10.2e ',epsrel);
    plot(t,fac*y(:,1),t,fac*y(:,3),'r',LW,1.0);
    set(gca,'Xlim',[0 tstopp]);
    grid ;
    xlabel('t(s)',FW,'Bold',FS,14);
    ylabel('\theta(degrees)',FW,'Bold',FS,14);
    title(st,FS,13);
    waitbar(k/antall);
end
close(h);

%==============================================================
function dydt = fcn(t,y)
% Equations of motion
global g m1 m2 l1 l2;
dydt = zeros(size(y));
dydt(1) = y(2);
dydt(2) = m2*l2*sin(y(3) - y(1))*y(4)^2 - g*(m1+m2)*sin(y(1));
dydt(3) = y(4);
dydt(4) = -l1*sin(y(3) - y(1))*y(2)^2 - g*sin(y(3));
%===============================================================
function M = mass(t,y)
% Mass-matrix 
global g m1 m2 l1 l2;
n = length(y);
M = zeros(n,n);
M(1,1) = 1.0;
M(2,2) = l1*(m1 + m2);
M(2,4) = m2*l2*cos(y(3) - y(1));
M(3,3) = 1.0;
M(4,2) = l1*cos(y(3) - y(1));
M(4,4) = l2;
%=====================================