% Program newton
% Computes the solution of Newton's
% 1. order equation (1671):
% dy/dx = 1-3*x + y + x^2 +x*y , y(0) = 0
% using Heun's method
%
clear all; close all;%clear all variables (and globals), delete all figures
clc; %clear command window

xend = 2;
steps=21;
x=linspace(0,xend,steps);
dx=x(2)-x(1)
y = zeros(steps,1); ya = y;  % allocate space
% y(1)= 0.0 ; x(1) = 0.0;
%fprintf('       t(s)       v(m/s)        Re \n\n');


for n = 1:steps - 1
   yp = y(n) + dx*f(x(n),y(n));                          % Predictor 
   y(n+1) = y(n) + 0.5*dx*(f(x(n),y(n)) + f(x(n+1),yp)); % Corrector
end


%% === Analytical solution
a = sqrt(2)/2;
t1 = exp(x.*(1+ x/2));
t2 = erf((1+x)*a)-erf(a);
ya = 3*sqrt(2*pi*exp(1))*t1.*t2 + 4*(1-t1)-x;

%% === Newton's solution
ynt = x - x.*x + x.^3/3 - x.^4/6; 

%% Plot results
h=plot(x,ya,x,y,'-.',x,ynt,'--');


%% Improve plot: set line width, labels, legends, fonts and title
set(h(:),'linewidth',2);

grid on

%% Labels and title
xlabel('x');ylabel('y','Rotation',0);
h3=legend('y_{a}','y_{heun}','y_{newton}'); set(h3,'box','off');
title('Newton''s equation (1671)');
