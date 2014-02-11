function val = fcng(alfa,k)
% Beregner funksjonen
% g = k - dn(alfa,k)/cn(alfa,k)/sqrt(2)
% der cn og dn er Jacobi-elliptiske funksjoner

tol = 1.0e-10; m = k^2;
[sn,cn,dn]= ellipj(alfa,m,tol);
val = k - dn/(sqrt(2)*cn);

