 function dpendel
% Computing the motion of a double pendulum 
% with unrestricted amplitudes.
% See Irgens : "Dynamikk", Chap. 9, examp.9.6
% Programmer: Jan B. Aarseth,
% Institute of Structural Engineering, NTNU

% ====== Variables ======
% The x-axis is horisontal with the positive
% direction to the right. The y-axis is vertical and
% positive in the downward direction.
% l1 is the length of the shoulder and m1 is the shoulder mass.
% l2 is the length of the elbow and m2 is the elbow mass.
% l1 makes the angle theta1 and l2 makes the angle theta2 with the y-axis.
% Both angles are positive in the counter-clock direction.
% 
% To write the equations as a system of 4 1.order equations,
% we put theta1 = y1 and theta2 = y3
% giving d(theta1)/dt = d(y1)/dt = y2
% and d(theta2)/dt = d(y3)/dt = y4
%
% Internal functions : FCN, MASS and PLOTTING 
%
global g m1 m2 l1 l2;
g = 9.81; % Gravity force

% ==== Data for masses and length of pendulums ===
% Standard-values: m1 = 0.5, m2 = 2.0, l1 = 1.5, l2 = 1.0
m1 = 0.5 ; m2 = 2.0;
l1 = 1.5; l2 = 1.0;

% === Reading time-data ====
tstop = input('Length of simulation(sec) = ? ');
tsteps = input('Number of time steps = ? ');
tstart = 0.0; 
tspan = linspace(tstart,tstop,tsteps);

% === Reading initial values of angles ====
theta1g = input('theta1(degrees) = ? ');
theta1 = theta1g*pi/180; % theta1 in radians
theta2g = input('theta2(degrees) = ? ');
theta2 = theta2g*pi/180; % theta2 in radians

% === Initial values of angular velocities ===
dtheta1 = 0.0; % d(theta)/dt for shoulder
dtheta2 = 0.0; % d(theta)/dt for elbow
L = l1 + l2;
FW = 'FontWeight'; FS = 'FontSize'; LW = 'LineWidth';

% ===== Plotting of initial values ====
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
st1 = sprintf('\\theta_1 = %3.1f\\circ ',theta1g);
xt = -L*0.9; yt = -L*0.9;
text(xt,yt,st1,FS,14,FW,'Bold');
st2 = sprintf('\\theta_2 = %3.1f\\circ ',theta2g);
yt = - L*0.7; 
text(xt,yt,st2,FS,14,FW,'Bold');
st = ('Press a key ....');
xt = - L*0.9; yt = L*0.9;
text(xt,yt,st,'Color',[0.5 0 1],FS,13);
pause;

% ==== Solving the differential equations ===
y0 =[theta1; dtheta1; theta2; dtheta2];
options = odeset('RelTol',1.0e-6,'Mass',@mass);
[t,y] = ode45(@fcn,tspan,y0,options);

% === Plotting the angles as function of time ===
clf reset;
fac = 180/pi;
plot(t,fac*y(:,1),t,fac*y(:,3),'r',LW,1.5);
grid ;
xlabel('t(s)',FW,'Bold',FS,14);
ylabel('\theta(degrees)',FW,'Bold',FS,14);
st = ('Blue : {\theta_1} , Red : {\theta_2} ');
title(st,FS,13);

pause;
close; clc;
bool1 = 1; bool2 = 1;
while (bool1)&(bool2)
   svar = input('CONTINUE WITH MORE GRAPHICS ?(y/n) ','s');
   svar = lower(svar);
   bool1 = ~strcmp(svar,'n');bool2 = ~strcmp(svar,'y');
end
if ~bool1
    return
end

% === Compute x and y-coordinates of the masses ===
n = length(t); n2 = 2*n;
X = zeros(n2,1); Y = zeros(n2,1);
X(1:2:n2) = l1*sin(y(1:n,1));
X(2:2:n2) = l1*sin(y(1:n,1)) + l2*sin(y(1:n,3));
Y(1:2:n2) = l1*cos(y(1:n,1));
Y(2:2:n2) = l1*cos(y(1:n,1)) + l2*cos(y(1:n,3));
X = [0;X]'; Y = [0;Y]';

% === Plotting motion of double pendulum ===
gstring1 = 'Xor'; gstring2 = 'Xor';
plotting(X,Y,L,gstring1,gstring2);
pause(3);
gstring1 = 'None'; gstring2 = 'None';
plotting(X,Y,L,gstring1,gstring2);
pause(3);
gstring1 = 'Xor'; gstring2 = 'None';
plotting(X,Y,L,gstring1,gstring2);

% ==============================================
function plotting(X,Y,L,gstring1,gstring2)
%=== Basic plotting routine ===
clf reset;
n = length(X) - 1;
axes('DrawMode','Normal','Box','on');
axis([-L L -L L],'ij');
set(gca,'DataAspectRatio',[1 1 1]);
lhandle1 = line(X(1:2),Y(1:2));
lhandle2 = line(X(1:2),Y(1:2));
set(lhandle1,'LineWidth',2,'EraseMode',gstring1);
set(lhandle2,'LineWidth',2,'EraseMode',gstring2,'Color',[1 0 0],'Marker','o');
for counter = 1: 2: n
    set(lhandle1,...
        'XData',[X(1) X(counter +  1)],...
        'YData',[Y(1) Y(counter +  1)]);
    set(lhandle2,...
        'XData',[X(counter + 1)  X(counter +  2)],...
        'YData',[Y(counter + 1)  Y(counter +  2)]);
    drawnow;
    pause(0.1)
end

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