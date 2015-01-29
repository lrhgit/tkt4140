function v = fskinit(eta)   
% Pohlhausen polynomials give initial values for the Falkner-Skan equation
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
% --------------------------------------------------------------------------
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