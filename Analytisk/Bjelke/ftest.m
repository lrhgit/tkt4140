% program ftest
% Beregner funksjonen
% g = k*sqrt(2) - dn(alfa,k)/cn(alfa,k)
% der cn og dn er Jacobi-elliptiske funksjoner
% fra ligningen k*sqrt(2) = dn(alfa,k)/cn(alfa,k)
%
alfa = 4.55336;
kvec = linspace(0.99966,0.99996,50);
n = length(kvec);
tol = 1.0e-6; 
for l = 1:n
    k = kvec(l);
    m = k^2;
    [sn,cn,dn]= ellipj(alfa,m,tol);
   val(l) = k*sqrt(2)- dn/cn;
end
plot(kvec,val);
shg

