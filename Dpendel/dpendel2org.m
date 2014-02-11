 function dpendel2org
% Løser bevegelsesligningene for en 
% dobbel-pendel med store utslag.
% Se Irgens : "Dynamikk", Kap. 9, eks.9.6
% Programmert av Jan B. Aarseth, NTNU
% Denne versjonen beregner pendelbevegelsen
% som funksjon av den relative nøyaktigheten RelTol
% i ligningsløseren.

% === Variable ===
% Setter theta1 = y1 og theta2 = y3
% slik at d(theta1)/dt = d(y1)/dt = y2
% og d(theta2)/dt = d(y3)/dt = y4
% Bruker interne funksjoner fcn, mass og plotting 
%
global g m1 m2 l1 l2;
g = 9.81; % Tyngden

% ==== Data for masser og pendel-lengder ===
m1 = 0.5 ; m2 = 2.0;
l1 = 1.5; l2 = 1.0;

% === Innlesning av tidsverdier ====
tstopp = 35;
antall = 350;
tstart = 0.0; 
tspan = linspace(tstart,tstopp,antall);

% === Innlesning av startvinkler ====
theta1g = 45;
theta1 = theta1g*pi/180; % theta1 i radianer
theta2g = 175;
theta2 = theta2g*pi/180; % theta2 i radianer

% === Andre startverdier ===
dtheta1 = 0.0; % Vinkelhastighet for l1
dtheta2 = 0.0; % Vinkelhastighet for l2
L = l1 + l2;
FW = 'FontWeight'; FS = 'FontSize'; LW = 'LineWidth';

% ===== Plotting av startverdier ====
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
text(-2,-1.75,st1,FS,14,FW,'Bold');
st2 = sprintf('\\theta_2 = %3.1f\\circ ',theta2g);
text(-2,-1.25,st2,FS,14,FW,'Bold');
pause;
% ==== Løser diff-ligningene ===
epsrel = 1.0e-6;
clf;
fac = 180/pi;
y0 =[theta1; dtheta1; theta2; dtheta2];
for k = 1:4
    epsrel = epsrel/10;
    options = odeset('RelTol',epsrel,'Mass',@mass);
    [t,y] = ode113(@fcn,tspan,y0,options);
    % === Plotter vinklene som funksjon av tiden ===
    st = sprintf('RelTol = %10.2e ',epsrel);
    figure(k);
    plot(t,fac*y(:,1),t,fac*y(:,3),'r',LW,1.5);
    grid ;
    xlabel('t(s)',FW,'Bold',FS,14);
    ylabel('\theta(grader)',FW,'Bold',FS,14);
    title(st,FS,13);
end
%==============================================================
function dydt = fcn(t,y)
% Differential-ligning
global g m1 m2 l1 l2;
dydt = zeros(size(y));
dydt(1) = y(2);
dydt(2) = m2*l2*sin(y(3) - y(1))*y(4)^2 - g*(m1+m2)*sin(y(1));
dydt(3) = y(4);
dydt(4) = -l1*sin(y(3) - y(1))*y(2)^2 - g*sin(y(3));
%===============================================================
function M = mass(t,y)
% Massematrisa 
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