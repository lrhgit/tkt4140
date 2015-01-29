clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

beta = -0.0;

[x,y]=falknerSkanSecant(beta);
xmax=max(x);

[x2,y2] = falknerSkanPenta(beta,0.1,xmax);

[x2,y2] = falknerSkanPenta(beta,h,xmax);


plot(y(:,2),x,'o',y(:,3),x,'-.',y2(1,1:end-1)',x2(1:end),'o',y2(2,1:end-1)',x2(1:end),'o');


ylabel('\eta')
xlabel('f'' , f"');
title('Solution of the Falkner-Skan equation')
legend('f'' Secant','f" Secant','f'' Penta','f" Penta')
