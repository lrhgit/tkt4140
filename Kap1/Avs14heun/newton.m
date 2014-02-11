% Program newton
% Computes the solution of Newton's
% 1. order equation (1671):
% dy/dx = 1-3*x + y + x^2 +x*y , y(0) = 0
% using Heun's method
%
clear all; % clear all variables (and globals),
close all; % delete all figures
clc;       % clear command window

xend = 2;
dx = 0.1;
steps = round(xend/dx) + 1;
y = zeros(steps,1); x = y ;ya = y ;ynt = y;  % allocate space
y(1)= 0.0 ; x(1) = 0.0;
%fprintf('       t(s)       v(m/s)        Re \n\n');

time1 = cputime;

for n = 1:steps - 1
   x(n+1) = n*dx;
   xn = x(n);
   fn = 1 + xn*(xn-3) + y(n)*(1 + xn);
   yp = y(n) + dx*fn; % Predictor
   xnp1 = x(n+1);
   fnp1 = 1 + xnp1*(xnp1-3) + yp*(1 + xnp1);
   y(n+1) = y(n) + 0.5*dx*(fn  + fnp1);
end

time2=cputime-time1;
fprintf('Elapsed time %f',time2);

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

hh(1,:)=xlabel('x');
hh(2,:)=ylabel('y','Rotation',0);
hh(3,:)=legend('y_{a}','y_{heun}','y_{newton}');

FS = 14;
set(hh(:),'FontName','Arial');
set(hh(:),'FontSize',FS);
set(gca,'FontSize',FS);
set(gca,'FontName','Arial');
set(hh(3),'box','off');
title('Newton''s equation (1671)')
%title('Newton''s equation (1671)','FontSize',FS)