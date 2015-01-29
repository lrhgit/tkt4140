% ================ Program btest ===========================
% This program computes the solution of the Blasius
% equation using finite differences and the Douglas algorithm:
%            f''' + f*f'' = 0
%         with boundary conditions :
%                f(0) = f'(0) = 0
%          f'(eta) -> 1.0, eta  -> infinity
% Input-variables : deta and iprt where deta is the grid-
% spacing and iprt is an output-control parameter.
% iprt = n , n > 0, prints results for every n gridpoint
%                   including the last point.    
% iprt = 0 prints output for only the first and last gridpoint
% ============================================================
clear
etainf = 5.8;
deta = 0.2;
nval = round(etainf/deta); % Number of intervals
etainf = nval*deta;
neq = nval - 1 ; % Number of equations
% ---- Allocation and initial values----
c2 = ones(neq,1); a2 = -c2; b1 = 2*deta*a2;
d1 = zeros(neq,1); a3 = d1; b4 = d1; c3 = d1; d2 = d1;
f = zeros(neq,1); g = f; u = f; v = f;
itmax = 15; epsi = 1.0e-5;

% === Initialize f and f'= g using Pohlhausen polynomials
t1 = deta/etainf;
for k = 1 : neq
   ksi = t1*k; ksi2 = ksi*ksi;
   f(k) = etainf*ksi2*(1.0 +ksi2*(ksi*0.2 -0.5));
   g(k) = ksi*(2.0 + ksi2*(ksi - 2.0));
end
fac = deta*0.5; 
for it = 1:1
   for k = 2 : neq-1;
       a3(k) = 1-fac*f(k);
       b4(k) = fac*(g(k+1) - g(k-1));
       c3(k) = 1 + fac*f(k);
       d2(k) = fac*f(k);
   end
   b4(1) = fac*g(2);
   c3(1) = 1 + fac*f(1);
   d2(1) = fac*f(1);
   b4(neq) = fac*(1 - g(neq-1));
   d1(neq) = -deta;
   d2(neq) = -1;
   [u,v] = douglas(a2,a3,b1,b4,c2,c3,d1,d2); % Solve system of equations
   g = u;
   f = v;
end % end of iteration loop




