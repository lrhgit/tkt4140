function tank
%   Løser vanntank-problemet
%  gitt i eksempel 2.6 i kompendiet.
%
global beta;
ny = 0.25 ; R = 2.75; H = 3.65;
t0 = 0.275;
beta = (3*(1-ny^2)*H^4/(R*t0)^2)^0.25;

solinit = bvpinit(linspace(0,1,10),@tankinit);
options = bvpset('RelTol',1.0e-3,'stats','on');
sol = bvp4c(@tankode,@tankbc,solinit,options);

xv = sol.x;
w  = sol.y;

%--- Plotting  av m(x)=-w'' og v(x) =-w'''
clf reset
plot(xv,-w(3,:)/beta,'k',xv,-w(4,:)/beta^2,'k-.');
%axis([0 1.4 0 etainf]);
title('Water tank')
ylabel('w')
xlabel('x')
grid
shg
fprintf('beta = %8.4f \n\n',beta);
xtab = [0:0.1:1];
tabell = deval(sol,xtab);
[xtab' tabell'];

% -----------------------------------------------

function dwdx = tankode(x,w)
global beta;
dwdx = zeros(4,1);
fac = -4*beta^4;
dwdx(1) = w(2);
dwdx(2) = w(3);
dwdx(3) = w(4);
dwdx(4) = fac*(w(1) + 1 - x);

% -----------------------------------------------
function res = tankbc(w0,w1)
% Boundary conditions : 
% w(0 = 0 , w'(0) = 0 , w''(1) = 0 w'''(1) = 0
res = [w0(1); w0(2); w1(3); w1(4)];
% -----------------------------------------------
function v = tankinit(x)   
global beta;
v = zeros(4,1);
fac = exp(-beta*x);
v(1) = fac*(sin(beta*x) + cos(beta*x));
v(2) = -2*beta*fac*sin(beta*x);
v(3) = 2*beta^2*fac*(sin(beta*x) - cos(beta*x));
v(4) = 4*beta^3*fac*cos(beta*x);

% ------------------------------------------------

