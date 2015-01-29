function tankv
%   Løser vanntank-problemet
%  gitt i eksempel 2.6 i kompendiet.
%  Veggtykkelsen varierer fra t = t0
%  ved x = 0 til t = t1 ved x = 1
%
global alpha beta;
ny = 0.25 ; R = 2.75; H = 3.65;
t0 = 0.275; t1 = 0.275;
beta = (3*(1-ny^2)*H^4/(R*t0)^2)^0.25;
alpha = (t0 - t1)/t0;

solinit = bvpinit(linspace(0,1,10),@tankinit);
options = bvpset('RelTol',1.0e-3,'stats','on');
sol = bvp4c(@tankode,@tankbc,solinit,options);

xv = sol.x;
w  = sol.y;

%--- Plotting  av m(x)=-w'' og v(x) =-w'''
% clf reset
% plot(xv,-w(3,:)/beta,'k',xv,-w(4,:)/beta^2,'k-.');
% %axis([0 1.4 0 etainf]);
% title('Water tank')
% ylabel('w')
% xlabel('x')
% grid
% shg
fprintf('beta = %8.4f \n\n',beta);
xtab = [0:0.05:1];
tabell = deval(sol,xtab);
[xtab' tabell']

% -----------------------------------------------

function dwdx = tankode(x,w)
global alpha beta;
dwdx = zeros(4,1);
fac = -4*beta^4;
fac1 = 1/(1 - alpha*x); fac2 = fac1^2;
dwdx(1) = w(2);
dwdx(2) = w(3);
dwdx(3) = w(4);
tmp1 = 6*alpha*fac2*(w(4) - fac1*alpha*w(3));
tmp2 = fac*fac2*(w(1) + (1-x)*fac1);
dwdx(4) = tmp1 + tmp2;

% -----------------------------------------------
function res = tankbc(w0,w1)
% Boundary conditions : 
% w(0 = 0 , w'(0) = 0 , w''(1) = 0 w'''(1) = 0
res = [w0(1); w0(2); w1(3); w1(4)];
% -----------------------------------------------
function v = tankinit(x)   
global alpha beta;
v = zeros(4,1);
fac = exp(-beta*x);
v(1) = fac*(sin(beta*x) + cos(beta*x));
v(2) = -2*beta*fac*sin(beta*x);
v(3) = 2*beta^2*fac*(sin(beta*x) - cos(beta*x));
v(4) = 4*beta^3*fac*cos(beta*x);

% ------------------------------------------------

