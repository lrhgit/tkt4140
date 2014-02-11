% =================== blaplot =========================
% Beregner og plotter fi som funksjon av s= f''(0) ved 
% bruk av skyteteknikk for l�sning av Blasius ligning
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20); %Default values for plotting.

sstart = 0.05; send = 0.8 ; antall = 30;
s = linspace(sstart,send,antall);
fi = s;
options= odeset('RelTol', 1.0e-5);
xspan = [0 5.75];
N=100;
dx=(max(xspan)-min(xspan));


for n = 1:antall
   y0 = [0.0 0.0 s(n)];
   [x,y] = ode45(@blasius,xspan,y0,options);
  % [t1,yrk] = rkn3('fcnbla',xspan,y0,dx,'rk4',0);
   fi(n) = y(end,2) - 1.0 ;
end
figure(1)
clf
h=plot(s,fi);
grid



%% Improve plot: set line width, labels, legends, fonts and title
grid on

xlabel('s');
ylabel('\phi');
title('Blasius - Nullpunkt for \phi')
