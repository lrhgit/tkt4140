 function dpendel
% Løser bevegelsesligningene for en 
% dobbel-pendel med store utslag.
% Se Irgens : "Dynamikk", Kap. 9, eks.9.6
% Programmert av Jan B. Aarseth, NTNU

% === Variable ===
% Setter theta1 = y1 og theta2 = y3
% slik at d(theta1)/dt = d(y1)/dt = y2
% og d(theta2)/dt = d(y3)/dt = y4
% Bruker interne funksjoner fcn, mass og plotting 
%
global g m1 m2 l1 l2;
g = 9.76; % Tyngden

% ==== Data for masser og pendel-lengder ===
% Standard-verdier: m1 = 0.5, m2 = 2.0, l1 = 1.5, l2 = 1.0
% m1 = 13680 ; m2 = 4560;
% l1 = 4.88; l2 = 4.88;
m1 = 0.5 ; m2 = 2.0;
l1 = 1.5; l2 = 1.0;

% === Innlesning av tidsverdier ====
tstopp = input('Simuleringstid = ? ');
antall = input('Antall tidskritt = ? ');
tstart = 0.0; 
tspan = linspace(tstart,tstopp,antall);

% === Innlesning av startvinkler ====
theta1g = input('theta1(grader) = ? ');
theta1 = theta1g*pi/180; % theta1 i radianer
theta2g = input('theta2(grader) = ? ');
theta2 = theta2g*pi/180; % theta2 i radianer

% === Andre startverdier ===
dtheta1 = 0.0; % Vinkelhastighet for l1
dtheta2 = 0.0; % Vinkelhastighet for l2
L = l1 + l2;

% ===== Plotting av startverdier ====
clf reset;
x1 = l1*sin(theta1); y1 = l1*cos(theta1);
x2 = x1 + l2*sin(theta2); y2 = y1 + l2*cos(theta2);
xl1 = [0; x1]; yl1 = [0 ; y1];
xl2 = [x1; x2]; yl2 = [y1 ; y2];
plot(xl1,yl1,xl2,yl2,'-ro','Linewidth',3);
axis([-L L -L L],'ij');
grid on
xlabel('x','FontWeight','Bold','FontSize',14);
ylabel('y','Rotation',0,'FontWeight','Bold','FontSize',14);
st = sprintf('l_1 = %3.2f , l_2 = %3.2f ,  m_1 = %3.2f ,  m_2 = %3.2f ',l1,l2,m1,m2);
title(st,'FontSize',13);
pause;

% ==== Løser diff-ligningene ===
y0 =[theta1; dtheta1; theta2; dtheta2];
options = odeset('RelTol',1.0e-4,'Mass',@mass);
[t,y] = ode45(@fcn,tspan,y0,options);

% === Beregner x og y-koordinater ===
n = length(t); n2 = 2*n;
X = zeros(n2,1); Y = zeros(n2,1);
X(1:2:n2) = l1*sin(y(1:n,1));
X(2:2:n2) = l1*sin(y(1:n,1)) + l2*sin(y(1:n,3));
Y(1:2:n2) = l1*cos(y(1:n,1));
Y(2:2:n2) = l1*cos(y(1:n,1)) + l2*cos(y(1:n,3));
X = [0;X]'; Y = [0;Y]';

% === Plotter ===
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
%=== Plotting ===
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