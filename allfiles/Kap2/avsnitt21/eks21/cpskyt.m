%===================== cpskyt ===============================
% The program solves the equation for the Couette-Poiseuille 
% flow using a shooting technique.
% The equation is given by:
%     u(0) = 0, u(1) = 1
%
%     u''(y) = - P
% Using rk4c to solve the system
%
% Calling function fcncp
%============================================================
clear all; close all; clc;
FS = 20; LW=3; set(0,'DefaultLineLineWidth',LW,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);

global P % Non-dimensional pressure gradient
dy = 0.05;
n = 1.0/dy; u = zeros(n+1,1); ua = u;y = u;
z = zeros(2,1);
P = input('Pressure gradient P = ? ');
% ==== Guessing initial values of s0 and s1:
s0 = input('s0 = ?');
s1 = input('s1 = ?');
fprintf('Pressure gradient P  = %10.3e \n',P);
fprintf('Initial values: s0 = %7.3f  s1 = %7.3f \n',s0,s1);

%% === Compute fi0 ===
z(1) = 0 ; z(2) = s0;
for k = 1: n
   y = (k-1)*dy;  
   z = rk4c(@fcncp,y,z,dy);
   %z = ode45(@fcncp,y,z,dy);
end
fi0 = z(1) - 1;
fprintf('  fio = %10.3e \n',fi0);
%%=== Compute fi1 ===
 z(1) = 0.0 ; z(2) = s1 ;
for k = 1: n 
   y = (k-1)*dy;  
   z = rk4c(@fcncp,y,z,dy);
end
fi1 = z(1) -1;
fprintf('  fi1 = %10.3e \n',fi1);
% === Compute s* ===
sstar = (fi1*s0 - s1*fi0)/(fi1 - fi0) ;
fprintf('  s*  = %10.3e \n',sstar);

%% === Final computation using s* ===
u(1) = 0; z(1) = 0.0 ; z(2) = sstar ;
for k = 1: n
   y = (k-1)*dy;  
   z = rk4c(@fcncp,y,z,dy);
   u(k+1) = z(1);
end
y = (0 : dy : 1.0)';
% --- Analytical solution ua 
ua = y.*(1 + P*(1-y)/2);
% === Table ===
% fprintf('\n    y           u        u-analyt.   \n\n'); 
% fprintf(' %7.3f   %10.4e  %10.4e \n',[y u ua]'); 

figure(1)
h=plot(y,u,'--',y,ua,'r:');



%% labels, legends and title
grid on

hh(1,:)=xlabel('y');
hh(2,:)=ylabel('u');
hh(3,:)=legend('u shooting','analytical');

% FS = 20;
% set(hh(:),'FontName','Arial');
% set(hh(:),'FontSize',FS);
% set(gca,'FontSize',FS);
set(hh(3),'box','off');
title('Shooting with RK4C')


