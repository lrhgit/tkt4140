% Program disserr
% Computes the dissipation error epsd of the 
% FTCS-scheme used on the diffusion(heat) equation.
% epsd = abs(G)/abs(Ga) where G is the numerical amplification 
% factor and Ga is the analytical amplification factor.
% r = dt/(dx)**2 (nondimensional form)
% G = 1 - 4*r*(sin(d/2))^2^, Ga = exp(-d^2*r), 0<= d <= pi)
% for stability : r <= 0.5

clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

hold on
axis([0 180 0 2])
for r = [1/6 0.25 0.5 ]
    dvec = linspace(0,pi)';
    n = length(dvec);
    epsd = zeros(n,1);
    for k = 1:n
        d = dvec(k);
        G = 1 - 4*r*sin(d*0.5)^2;
        Ga = exp(-d^2*r);
        epsd(k) = abs(G/Ga);
    end
    plot(dvec*180/pi,epsd,'k')
end
hold off
grid on
shg