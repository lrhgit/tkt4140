% =================== blaplot =========================
% Beregner og plotter fi som funksjon av s= f''(0) ved 
% bruk av skyteteknikk for l�sning av Blasius ligning
clear all; close all; clc;

sstart = 0.05; send = 0.8 ; antall = 30;
s = linspace(sstart,send,antall);
fi = s;
options= odeset('RelTol', 1.0e-5);
xspan = [0 5.75];
for n = 1:antall
   y0 = [0.0 0.0 s(n)];
   [x,y] = ode45(@blasius,xspan,y0,options);
   fi(n) = y(end,2) - 1.0 ;
end
figure(1)

h=plot(s,fi);
grid



%% Improve plot: set line width, labels, legends, fonts and title
set(h(:),'linewidth',2);

grid on

hh(1,:)=xlabel('s');
hh(2,:)=ylabel('\phi');

FS = 20;
set(hh(:),'FontName','Arial');
set(hh(:),'FontSize',FS);
set(gca,'FontSize',FS);
%set(hh(3),'box','off');
title('Blasius - Nullpunkt for \phi','Fontsize',FS,'FontWeight','Bold')
