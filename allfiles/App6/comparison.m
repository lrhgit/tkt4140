clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

beta = -0.0;

[x,y]=falknerSkanSecant(beta);
xmax=max(x);

h=x(2)-x(1) % Compute h used for Secant method
disp('h - '), h

h=h/10.0


[x2,y2] = falknerSkanPenta(beta,h,xmax);


%plot(y(:,2),x,'o',y(:,3),x,'-.',y2(1,1:end-1)',x2(1:end),'o',y2(2,1:end-1)',x2(1:end),'o');
plot(y(:,2),x,y(:,3),x,y2(1,1:end-1)',x2(1:end),y2(2,1:end-1)',x2(1:end));


ylabel('\eta')
xlabel('f'' , f"');
title('Solution of the Falkner-Skan equation')
legend('f'' Secant','f" Secant','f'' Penta','f" Penta')
