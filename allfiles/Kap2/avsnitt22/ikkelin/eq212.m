% ============= Program eq212 ==============
%  Løser eksempel i avsnitt 2.2 i kompendiet.
%  
%  Ligningen er gitt ved:
%  y''(x) = a*y^2 , a > 0
%  y(0) = 6/a, y(1) = 3/(2*a)
%  I programmet er 1<= a <= 10
%  Finner løsningen yII ved skyteteknikk.
%  Den analytiske løsningen er gitt ved :
%  y(x) = (6/a)*P(x + C1,0,C2)
%  der P er Weierstrass elliptiske P-funksjon. 
% 
%  I programmet setter vi y(1)= y, y(2) = y'
% 
% s = y'(0) = y(2)(0)
% Bruker funksjonen fcn21.m
%
% Skriv inn en a-verdi mellom 1 og 10 
%
clear all; close all; clc; 

global a;
a = 2.0;
if (a > 10) || ( a < 1)
    error('Parameter a is out of range!');
end
xspan = [0 1];
y0 = 6/a;
y1 = 3/(2*a);
%% Beregn startverdier
s0 = -53.788/a;
s1 = s0*0.95;

%% Beregner fi0
ys = [y0 s0];
options = odeset('RelTol',1.0e-7);
[x,y] = ode45(@fcn21,xspan,ys,options);
fi0 = y(end,1) - y1 ;

%% Startverdier for iterasjonen
itmax = 10; epsi = 1.0e-5; it = 0; ds = 1;


%% Skriv overskrift for tabell
fprintf('        itr.      s          ds\n\n');

%%   Starter iterasjon
tic
while(abs(ds) > epsi) && (it < itmax)
   it = it + 1;
   ys = [y0 s1];
   [x,y] = ode45(@fcn21,xspan,ys,options);
   fi1 = y(end,1) - y1;
   ds = -fi1*(s1 - s0)/(fi1 - fi0);
   s = s1 + ds;
   s0 = s1;
   s1 = s;
   fi0 = fi1;
   fprintf('%10d %12.6f %12.3e\n',it,s,ds);
end
toc

%% --- Beregner en tabell for y og y' 
ys = [y0 s];
sol = ode45(@fcn21,xspan,ys,options);
xval = (0 : 0.05 :1);
tabell = deval(sol,xval);
fprintf('\n         x         y          y'' \n\n');
fprintf(' %12.2f %10.5f %10.5f \n',[xval' tabell']');

%% --- Plotter y
FS = 'FontSize'; FW = 'FontWeight';
plot(sol.x,sol.y(1,:)')
grid 
xlabel('x',FS,14,FW,'Bold')
ylabel('y','Rotation',0,FS,14)
st = sprintf('Eksempel i avsnitt 2.2 med a = %5.2f \n',a);
title(st,FS,13)
