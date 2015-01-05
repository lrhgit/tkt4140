% program analyticv2
% Som analyticv1, men bare for r = 0 -> w(0,t):
% Beregner den analytiske løsningen av diff. ligningen
%   dw/dt = w''(r) + w'(r)/r , w = w(r,t), 0 < r < 1

t = 0.2952; % Tid;
% Beregner 1. ledd forseg selv
n = 1;
lam1 = j0zero(n);
term1 = exp(-t*lam1^2)/(besselj(1,lam1)*lam1^3);  
sum1 = term1; sumtol = 1.0e-6; test = 1;
while test > sumtol
    n = n + 1;
    lamn = j0zero(n);
    term = exp(-t*lamn^2)/(besselj(1,lamn)*lamn^3);       
    sum1 = sum1 + term;
    test = abs(term/term1);
end
wzero = 8*sum1;
u = 1 - wzero