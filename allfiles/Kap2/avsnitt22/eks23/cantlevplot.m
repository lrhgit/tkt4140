% ======== Program cantlevplot ======== 
% Computes initial values for shooting
% technique in program cantilev1 and 2 
% Computes dtheta(1)/ds as a function of S where
% S = dtheta(0)/ds
% Here theta is the slope-angle and s the arc-length.
%  
clear all; close all; clc;
FS = 20; LW=3; set(0,'DefaultLineLineWidth',LW,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);


global alpha2 ; %% used by the beam function.

alpha2 = input('alpha2 = ? ');

Sstart = .1; Send =  3.5 ; antall = 20;
S = linspace(Sstart,Send,antall);
phi = S;

options= odeset('RelTol', 1.0e-5);
sspan = [0 1.0];

for n = 1:antall
   y0 = [0.0 S(n)];
   [s,y] = ode45(@beam,sspan,y0,options);
   phi(n)=y(end,2) ;
end


plot(S,phi)
grid
xlabel('s')
ylabel('\phi')
st = sprintf('Root of \\phi   \\alpha^2 = %4.2f',alpha2);
title(st)
