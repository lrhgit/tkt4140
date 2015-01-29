% Program gfaktor
% Computes G and Ga for the
% FTCS-scheme used on the diffusion(heat) equation.
% G is the numerical amplification factor
% and Ga is the analytical amplification factor.
% r = dt/(dx)**2 (nondimensional form)
% G = 1 - 4*r*(sin(d/2))^2^, Ga = exp(-d^2*r), 0<= d <= pi)
% for stability : r <= 0.5

clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);


hold all;
axis([0 180 0 1])
linespec=['k-.','k','b-.','b','r-.','r'];
fac = 180/pi;
i=1;
for r = [ 1/6  0.25  0.5 ]
   
    dvec = linspace(0,pi)';
    n = length(dvec);
    G = zeros(n,1);
    Ga = zeros(n,1);
    for k = 1:n
        d = dvec(k);
        G(k) = 1 - 4*r*sin(d*0.5)^2;
        Ga(k) = exp(-d^2*r);
    end
    %plot(dvec*fac,G,'k-.',dvec*fac,Ga,'k')
   % plot(dvec*fac,G,linespec(i:i+2),dvec*fac,Ga,linespec(i+3),'DisplayName',{strcat('G =',num2str(r)),strcat('Ga =',num2str(r))})
    plot(dvec*fac,G,linespec(i:i+2),dvec*fac,Ga,linespec(i+3),'DisplayName',{strcat('G =',num2str(r)),'G_a'});
    legend('-DynamicLegend');
    i=i+4; 
end
hold off
grid on
