% Program ribbe22
% Beregner eksemplet i avsnitt 3.2.1,
% tilfelle II med foroverdifferanser,
% med bruk av tdma.
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);


beta = 2;
h = 0.1;
n = round(1.0/h); % Antall ligninger
a = ones(n,1); % Underdiagonal
bfac = (beta*h)^2;
b = -(2 + bfac)*a; % Hovediagonal
b(1) = 2;
c = a; % Overdiagonal
c(1) = -(2.0 - bfac); 
d = zeros(n,1); % Hyre side
d(n) = -1;
thet = tdma(a,b,c,d);
thet = [thet; 1];
% analytisk lsning
xa = (0 : h: 1.0)';
theta = cosh(beta*xa)/cosh(beta);
fprintf('  x      theta     analyt \n');
fprintf('%4.2f  %9.5f  %9.5f \n',[xa thet theta]');

h=plot(xa,theta,xa,thet,':');

