%================ Program beamcola =========================
%                   Example 2.2A
% The program computes the deflection of a beam-column
% with a uniform lateral load of intensity q and an axial
% compressive force P.
%   Differential equation:
%   u''(x) + theta^2*u(x) = - theta^2*(1-x^2)/2, -1< x <1
%   u(-1) = 0 , u(1) = 0 
%   with theta^2 = P*L^2/EI. Buckling for theta = Pi/2
%
%  Analytical solution:
%  ua = (cos(theta*x)/cos(theta)-1)/theta2 -(1-x.^2)/2;
%
% The differential equation is solved by a shooting
% technique and integrated using RK4C
%
% Using function fcnbeama
%
%==========================================================
clear; clear global theta2 velg;
global theta2 velg

% === Initialize ===
dx = 0.1;
n = 2.0/dx; u = zeros(n+1,1); u0 = u; u1 = u; ua = u;
y = zeros(2,1);
theta = 1.0   ;  % Load-parameter
% Note : Buckling for theta = Pi/2
theta2 = theta^2;
fprintf('  Load-parameter theta  = %10.3e \n',theta);

% === Compute u0 ===
velg = 1; % Switch to first equation in function fcnknekk
u0(1)= 0.0; y(1) = 0.0 ; y(2) = 0.0;
for k = 1: n
   x = (k-1)*dx - 1;  
   y = rk4c(@fcnbeama,x,y,dx);
   u0(k+1) = y(1);
end

% === Compute u1 ===
velg = 2; % Switch to second equation in function fcnbeama
u1(1)= 0.0; y(1) = 0.0 ; y(2) = 1.0 ;
for k = 1: n 
   x = (k-1)*dx - 1;  
   y = rk4c(@fcnbeama,x,y,dx);
   u1(k+1) = y(1);
end

% === Complete solution u ===
u = u0 - (u0(n+1)/u1(n+1))*u1;
x = (-1 : dx : 1.0)';
% --- Analytical solution ua 
ua = (cos(theta*x)/cos(theta)-1)/theta2 -(1-x.^2)/2;

% === Output ===
fprintf('\n    x       u-comput.    u-analyt.   \n\n'); 
fprintf(' %7.3f   %10.4e  %10.4e \n',[x u ua]'); 


