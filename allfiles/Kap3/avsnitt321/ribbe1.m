% Program ribbe1
% Beregner eksemplet i avsnitt 3.2.1
% , tilfelle I, med bruk av tdma.
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

beta = 5;
h = 0.1;
n = round(1.0/h) - 1; % Antall ligninger
a = ones(n,1); % Underdiagonal
b = -(2 + (beta*h)^2)*a; % Hovediagonal
c = a; % Overdiagonal
d = zeros(n,1); % Hyre side
d(n) = -1;

thet = tdma(a,b,c,d);
thet = [0;thet; 1];

% analytisk lsning
xa = (0 : h: 1.0)';
theta = sinh(beta*xa)/sinh(beta);
fprintf('  x      theta     analyt \n');
fprintf('%4.2f  %9.5f  %9.5f \n',[xa thet theta]');

h=plot(xa,theta,xa,thet,':');
