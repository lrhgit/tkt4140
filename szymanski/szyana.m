% program szyana
clear
s = 0;
r = 0.0; t = 0.01;
test = 1.0; tol = 1.0e-7;
k = 0;
while (test > tol)
    k = k + 1;
    lk = j0zero(k);
    arg = lk*r;
    term = besselj(0,arg)/(besselj(1,lk)*lk^3);
    fac = exp(-lk^2*t);
    s = s + term*fac;
    test = abs(term);
end
w = 8*s
u = 1 - r^2 - w