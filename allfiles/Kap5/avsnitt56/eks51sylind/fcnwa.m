function w = fcnwa(r,t)
% Beregner den analytiske løsningen av diff. ligningen
%   dw/dt = w''(r) + w'(r)/r , w = w(r,t), 0 <= r <= 1
% der t er tiden.
% Nøyaktigheten kan kontrolleres med å forandre verdien
% av sumtol som angir den relative nøyaktigheten.
% Langsom konvergens ledd nær r = 1 for liten t;
%
if r == 1
    w = 0;
    return
end
if (r < eps) & (t < eps )
    w = 1;
    return
end
% === Beregner 1. ledd separat
n = 1;
lam1 = j0zero(n);
arg1 = r*lam1;
term1 = exp(-t*lam1^2)*besselj(0,arg1)/(besselj(1,lam1)*lam1^3);  
sum1 = term1; sumtol = 1.0e-8; test = 1;
while test > sumtol
    n = n + 1;
    lamn = j0zero(n);
    arg = r*lamn;
    term = exp(-t*lamn^2)*besselj(0,arg)/(besselj(1,lamn)*lamn^3);       
    sum1 = sum1 + term;
    test = abs(term/term1);
end
w = 8*sum1;
