% === Program eks22 ===
% Programmet løser randverdi-problemet i 
% avsnitt 2.2 i kompendiet ved bruk av skyteteknikk.
% Bruker sekantmetoden for nullpunktbestemmelse.
% Ligning : y''(x)= (3/2)*y(x)^2
%           y(0) = 4, y(1) = 1
% 
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20); %Default values for plotting.

%xspan = [0.0 1.0];
xmax=1.0;
xmin=0.0;
xspan = [xmin,xmax];

% ----Tipper to verdier s0 og s1 for y'(x0)
s0 = -35; s1 = -40; %% close to yII
%s0 = -3; s1 = -6;  %% close to yI

% ---- Beregner fi0
y0 = [4.0; s0]; % Startverdier
[x,y] = ode45(@fcn22,xspan,y0);
fi0 = y(end,1) -1;

% ----- Startverdier for iterasjonen
itmax = 10; epsi = 1.0e-4; it = 0; ds = 1;
options = odeset('RelTol',1.0e-1);
fprintf('        itr.      s          ds\n\n');

%% ===== Starter iterasjon =====
tic
while(abs(ds) > epsi) & (it < itmax)
   
   it = it + 1;
   y0 = [4.0; s1];
   [x,y] = ode45(@fcn22,xspan,y0,options);
   fi1 = y(end,1) -1;
   ds = -fi1*(s1 - s0)/(fi1 - fi0);
   s0 = s1;
   s1 = s1 + ds;
   fi0 = fi1;
   if abs(s1) > 1 
      ds = ds/s1;
   end;
   fprintf('%10d %12.6f %12.3e\n',it,s1,ds);
end
toc


%% Analytical solution
ya= 4./(1+x).^2;


%% Plot the results
% [ax,h1,h2]=plotyy(x,y(:,1),x,y(:,2));
% title('dy^2/x^2 = 3/2 y^2')
% ylabel('y');
% axes(ax(2)); ylabel('y'''); %Make ax(2) the current axes and set label

plot(x,y(:,1),x,ya,':');
title('dy^2/x^2 = 3/2 y^2')
ylabel('y');
legend('y','y_a');

%axes(ax(2)); ylabel('y'''); %Make ax(2) the current axes and set label


% %% ---- Skriv ut en tabell for y og y'
% fprintf('\n           x         y           y''\n\n');
% fprintf(' %12.2f %12.6f %12.6f\n',[x y]');
