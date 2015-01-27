function blasius
%   Solves the Blasius equation
%  'etainf' is a variable to facilitate experimentation.
%
global etainf;
etainf = 6;

solinit = bvpinit(linspace(0,etainf,10),@blainit);

options = bvpset('stats','on');

sol = bvp4c(@blaode,@blabc,solinit,options);
eta = sol.x;
f = sol.y;

fprintf('\n');
fprintf(' f''''(0) = %7.5f.\n',f(3,1))

clf reset
plot(f(2,:),eta);
axis([0 1.4 0 etainf]);
title('Blasius equation')
ylabel('\eta')
xlabel('df/d\eta')
shg

% --------------------------------------------------------------------------

function dfdeta = blaode(eta,f)
%   
dfdeta = [ f(2); f(3); -f(1)*f(3)];

% -------------------------------------------------------------------------
function res = blabc(f0,finf)
% Boundary conditions for Blasius equation
res = [f0(1); f0(2); finf(2) - 1];
% -------------------------------------------------------------------------
function v = blainit(eta)
%   
global etainf;
%  The guess for the solution satisfies the boundary conditions.
x = eta/etainf; x2 = x*x;
v = zeros(3,1);
v(1) = etainf*x2*(1 + x2*(0.2*x - 0.5));
v(2) = x*(2.0 + x2*(x - 2));
v(3) = (2 + x2*(4*x - 6))/etainf;

% --------------------------------------------------------------------------

