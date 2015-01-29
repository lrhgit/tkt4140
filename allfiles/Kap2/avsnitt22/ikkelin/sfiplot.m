% Beregner fi som funksjon av s % 
% for problemet i avsnitt 2.2.
%
% Ved å utvide s-området kan man se at finnes er to løsninger.
%
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

smax=-3.0; smin = -40.0 ; antall = 200;
s=linspace(smin,smax,antall);
fi = s; %intialize fi with same dim as s

options= odeset('RelTol', 1.0e-5);
xspan = [0 1.0];

for n=1:antall
   y0 = [4.0 s(n)];
   [x,y] = ode45('fcn22',xspan,y0,options);
   fi(n)=y(end,1) - 1.0 ;
end
figure(1)


%% Find zeros of fi
smallValue=0.05;
index=find(abs(fi)<smallValue);


%% Plot the results
plot(s,fi,s(index),fi(index),'ro');


xlabel('s');
ylabel('\phi');

grid
title('Nullpunkt for \phi','FontWeight','Bold')