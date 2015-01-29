%===================== vibrana ===================
% Large amplitude vibration of a mass on a frictonless
% foundation. The program computes the analytical solution
% of the equation: u''(t) + u(t)^3 = 0 with u(0) = -u0
% We have used u0 = 1.
% The complete period T = 7.416298. Only 1/2 of the period
% is computed because of the symmetry.
%======================================================
clear
u0 = 0.1;
t = [0 : 0.25 : 20.0]';
T = t*u0;
cost = -u0*cos(T);
[sn,cn,dn] = ellipj(T,0.5);
u = - cn*0.1;
fprintf('       t         u  \n\n');
fprintf(' %10.3f %15.4e \n',[t,u]');
clf
plot(t,u,t,cost,'-.')
xlabel('t','Fontsize',14,'Fontweight','Bold')
ylabel('u','Fontsize',14)
title('Large amplitude vibration','Fontsize',14)
legend('u','cost',2)   


