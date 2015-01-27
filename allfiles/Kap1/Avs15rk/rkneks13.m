% program rkneks13
clear all; close all; clc;
FS = 20; set(0,'DefaultLineLineWidth',3,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);
tspan = [0 10];
y0 = [0 ; 0];

dt = 1.;
g = 9.81; a = 0.007;

[t,ye]  = rkn('odefun',tspan,y0,dt,'rk1',g,a);
[t,yh]  = rkn('odefun',tspan,y0,dt,'rk2',g,a);
[t,yrk] = rkn('odefun',tspan,y0,dt,'rk4',g,a);

%% Plot the solution
h=plot(t,ye(:,2),t,yh(:,2),':',t,yrk(:,2));
grid on

%% Title, labels and legends
title('Falling sphere with constant C_{D}');
xlabel('t');  ylabel('z');
h3=legend('rk1','rk2','rk4'); set(h3,'box','off');

