function bladiff2
%  Solves the Blasius equation
%  for a given value of beta using the Matlab-function
%  bvp4c. 
%  This version uses nested functions.

etainf = 5.75;
solinit = bvpinit(linspace(0,etainf,10),@blainit);

% set 'stats' = 'on' to get information
options = bvpset('stats','off');

sol = bvp4c(@blaode,@blabc,solinit,options);
eta = sol.x;
f = sol.y;

FS = 'FontSize';
clf 
plot(f(2,:),eta,'k',f(3,:),eta,'k-.');
axis([0 1 0 etainf]);
title('Blasius equation',FS,14);
ylabel('\eta',FS,14,'Rotation',0);
xlabel('f''(\eta), f''''(\eta) \rightarrow ',FS,14)
grid
shg

% --------------------------------------------------------------------------
function dfdeta = blaode(eta,f)
% The Blasius equation
dfdeta = [ f(2); f(3); -f(1)*f(3)];
end
% -------------------------------------------------------------------------
function res = blabc(f0,finf)
% Boundary conditions for the Blasius equation
res = [f0(1); f0(2); finf(2) - 1];
end
% -------------------------------------------------------------------------
function v = blainit(eta)   
% Pohlhausen polynomials give initial values for the Blasius equation
% These polynomials satisfy the boundary conditions
% f(0) = 0 , f'(0)=0, f'(etainf) = 1
x = eta/etainf; x2 = x*x;
v = zeros(3,1);
v(1) = etainf*x2*(1 + x2*(0.2*x - 0.5));
v(2) = x*(2.0 + x2*(x - 2));
v(3) = (2 + x2*(4*x - 6))/etainf;
end
end
