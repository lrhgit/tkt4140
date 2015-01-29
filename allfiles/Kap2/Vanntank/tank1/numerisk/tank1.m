% Program tank1
% Section 2.4.2 : Cylindrical water-tank
% with a constant wall-thickness.
% Compute displacement, shearforce and bending moment using
% a shooting technique.
% Equation :
%  w''''(x) + 4*beta^4*w(x)= -4*beta^4*(1-x)
% Boundary conditions : w(0) = 0, w'(0) = 0, w''(1)=0, w'''(1) = 0
% The equation is linear but stiff for large values of beta.
%
clear; clear global beta4;
global beta4;
beta = input('beta = ? ');
fprintf('     beta = %7.3f\n',beta);
beta4 = beta^4;
xspan = (0: 1.0);
s = [0 0 1];  r = [0 1 0];
options = odeset('Reltol',1.0e-7,'AbsTol', 1.0e-20);
phi = zeros(3,1); psi = phi;
% ===== Skyter tre ganger for å finne s* og r* =====
for k = 1:3   
   y0 = [0.0; 0.0 ; s(k) ; r(k)]; 
   [x,y] = ode45(@fcntank1,xspan,y0,options);
   phi(k) = y(end,3);
   psi(k) = y(end,4);
end
nev = (psi(3)-psi(1))*(phi(2)-phi(1)) - (phi(3) - phi(1))*(psi(2)-psi(1));
rstar= (phi(3)*psi(1) - psi(3)*phi(1))/nev;
sstar = (psi(2)*phi(1) - phi(2)*psi(1))/nev;
fprintf('     r* = %12.5e  s* = %12.5e \n\n',rstar,sstar);
% ===== Table of displacement w,slope dw/dx etc. =====
xspan = (0: 0.1 :1.0);
y0 = [0.0 ;0.0 ;sstar; rstar];
[x,y] = ode45(@fcntank1,xspan,y0,options);
fprintf('       x         w             dw/dx           m(x)            v(x)\n\n');
fprintf( '%10.3f  %13.5e  %13.5e   %13.5e  %13.5e \n',[x y]');
m = -y(:,3)/beta; v = -y(:,4)/beta^2;
% ====== Plotting m(x)/beta and v(x)/beta^2 =====
clf
plot(x,m,'k-',x,v,'k-.','LineWidth',1.25);
grid on
xlabel('x','FontSize',14)
st = sprintf('Water tank. \\beta = %4.2f',beta);
title(st,'Fontsize',14)
legend('m(x)/\beta','v(x)/\beta^2')




  
 