%********************** Program heatex ***************************
%                                                                *
%   This program computes the steady state temperature           *
%   in a two-tube heat exchanger working in counterflow.         *
%                                                                *
%   The differential equations is given by :                     *
%                                                                *
%                 u(x)' = - aa*[u(x) - v(x)]                     *      
%                 v(x)' = - ab*[u(x) - v(x)]                     *
%                                                                *
%   where u is the temperature in the outer tube                 *
%         v is the temperature in the inner tube                 *
%                aa = h*p*l/(ma*ca)                              *
%                ab = h*p*l/(mb*cb)                              *
%                                                                *
%   Index a means the outer tube and index b  the inner tube.    *
%                                                                *
%     h = total heat transfer coefficient (W/m*m/deg. C)         *
%     p = perimeter of tube b (m)                                *
%     l = length of tubes (m)                                    *
%     ma = rate of mass in tube a (kg/s)                         *
%     ca = specific heat of medium in tube a. (J/kg/deg. C)      *
%                                                                *
%   Boundary conditions:                                         *
%                u(0) = 0,  v(1) = 1                             *
%                                                                *                                                                
%*****************************************************************
clear; close;
dx = 0.10;
if mod(1,dx) ~= 0
    error('No. of equations not an integer!')
end
x = (0 : dx: 1)';
n = 1/dx; % No. of equations
% --- Allocate space ---
a1 = zeros(n,1);
a3 = a1; b1 = a1; b2 = a1; b3 = a1;
b4 = a1; c2 = a1; c4 = a1; d1 = a1; d2 = a1; u = a1; v = a1;
aa = 0.9;
ab = 0.2;

fprintf('Grid-step dx = %8.3e \n\n',dx);
%
%--- Compute coefficients ----
%
for k = 1:n
    a1(k) =  dx*aa - 2;
    a3(k) =  dx*ab;
    b1(k) =  2 + dx*aa;
    b2(k) =  - dx*aa;
    b3(k) =    dx*ab;
    b4(k) =  - (dx*ab + 2);
    c2(k) =  - dx*aa;
    c4(k) =  2 - dx*ab;
end
d1(n) = dx*aa;
d2(n) = dx*ab -2;

%--- Compute u and v ----
      
[u,v]= bitris(a1,a3,b1,b2,b3,b4,c2,c4,d1,d2);

u = [0 ; u];
v = [v ; 1];
      
% fprintf('%6.2f %12.5e %12.5e \n',[x u v]');
% fprintf('\n');
%
%     ==== Compute the "analytical" solution ====
%
v0 = 0.87425359;
ua = 1 - exp(-0.7*x);
va = v0*(1 + 2*ua/7);
ua = v0*9*ua/7;
% fprintf('%6.2f %12.5e %12.5e \n',[x ua va]');
% fprintf('\n');
%
% --- Output in centigrades ---
% 
u = 100 - 70*u;
v = 100 - 70*v;
%fprintf('%6.2f %12.5e %12.5e \n',[x u v]');
clf
FS = 'FontSize';FW = 'FontWeight';
plot(x,u,'k',x,v,'k-.');
ylim([0 100]);
xlabel('x',FS,14,FW,'Bold')
ylabel('Temperature (\circC)',FS,14)  
title('Heat exhanger - counterflow',FS,14)
legend('lubricate','water')
     