function tval = txy(x,y)
% Computes the analytical solution of the Laplace  
% equation  Txx + Tyy = 0 on the unit square in figure 7.6    

s = 0.0; sn = 1;
ph = pi*0.5;
fac = cosh(ph*y)/(ph*cosh(ph));
n = 1; test = 1;
reltol = 1.0e-8;
while (test > reltol)
    lam = (2*n-1)*ph;
    term1 = cosh(lam*y)/(lam*cosh(lam));
    test = abs(term1/fac);
    ds = sn*term1*cos(lam*x);
    s = s + ds;
    n = n + 1;
    sn = -sn;
end
tval = 2*s;
