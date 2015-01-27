function test2
%  Solves the Falkner-Skan equation
%  for a given value of beta using Matlab-function
%  bvp5c. 
%  This version uses nested functions.
betasep = -0.19883768;
beta = input(' beta = ?');
while(beta < betasep) | (beta > 1.999)
   fprintf('beta = %7.3e is invalid!\n',beta);
   disp(' Try again !');
   beta = input(' beta = ?');
end
etainf = svalue(beta);
solinit = bvpinit(linspace(0,etainf,10),@fskinit);
% set 'stats' = 'on' to get information
options = bvpset('stats','on');
sol = bvp5c(@fskode,@fskbc,solinit,options);
eta = sol.x;
f = sol.y;
fprintf('\n');
fprintf(' f''''(0) = %7.5f.\n',f(3,1))
% Limits for plotting
v1 = ceil(10*f(3,1))*0.1;
v1 = max(1.1,v1);
v2 = ceil(10*etainf)*0.1;
FS = 'FontSize';
clf 
plot(f(2,:),eta,'k',f(3,:),eta,'k-.');
axis([0 v1 0 v2]);
stitle = sprintf('Falkner-Skan equation. \\beta = %5.3f',beta);
title(stitle,FS,14);
ylabel('\eta',FS,14,'Rotation',0);
xlabel('f''(\eta), f''''(\eta) \rightarrow ',FS,14)
grid
shg
% -------------------------------------------------------------------
function dfdeta = fskode(eta,f)
% The Falkner-Skan equation
dfdeta = [ f(2); f(3); -f(1)*f(3) - beta*(1 - f(2)^2)];
end
% -------------------------------------------------------------------
function res = fskbc(f0,finf)
% Boundary conditions for the Falkner-Skan equation
res = [f0(1); f0(2); finf(2) - 1];
end
% -------------------------------------------------------------------
function v = fskinit(eta)   
% Pohlhausen polynomials give initial values for the FSK-equation
% These polynomials satisfy the boundary conditions
% f(0) = 0 , f'(0)=0, f'(etainf) = 1
x = eta/etainf; x2 = x*x;
v = zeros(3,1);
fac = beta*etainf^2/6;
fac1 = fac*(0.5 - x + x2*(0.75 -0.2*x));
fac2 = fac*(1 - 3*x +x2*(3 - x ));
fac3 = fac*(1 + x2*(9-4*x)-6*x);
v(1) = etainf*x2*(1 + x2*(0.2*x - 0.5)+ fac1);
v(2) = x*(2.0 + x2*(x - 2)+ fac2);
v(3) = (2 + x2*(4*x - 6)+ fac3)/etainf;
end
% -------------------------------------------------------------------
function etamax = svalue(beta)
% Estimate of etamax for a given value of beta
% See appendix 3,part 3.
p1 = [0.633 -1.68 5.76] ; p2 = [36.76  2.0 5.87];
p3 = [0.125  -0.9 5.463]; 
if ( beta <= 1.0)
   if ( beta >= 0.0 ) % 0 <= beta <= 1.0
      p = polyval(p1,beta);
      etamax = round(10.0*p)*0.1;
   else
      p = polyval(p2,beta); % betasep <= beta < 0
      etamax = round(10.0*p)*0.1;
   end
else
   p = polyval(p3,beta); % 1 < beta <= 1.999
   etamax = round(10.0*p)*0.1;
end
end
end
