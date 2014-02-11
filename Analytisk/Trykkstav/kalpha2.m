function k = kalpha2(alpha)
% Elastica - Large deflection of a beam
% Computes k for a given value of alpha
% from the equation sn(alpha,k) = 1
% Range : pi/2 <= alpha <= 6.12778
 persistent kvec alphavec
if isempty(kvec) 
   kvec = [0.0 0.08716 0.17365 0.25882 0.34202 0.42262 0.5 0.57358 ...
           0.64279 0.70711 0.76604 0.81915 0.86603 0.90631 0.93969 ...
           0.96593 0.98481 0.99619 0.99905 0.99939 0.99966 0.99985 0.99996]'; 
   alphavec = [1.57080 1.57379 1.58284 1.59814 1.62003 1.64900 1.68575 1.73125 ...
          1.78677 1.85407 1.93558 2.03472 2.15652 2.30879 ...         
          2.50455 2.76806 3.15339 3.83174 4.52022 4.74272 5.02986 5.43491 6.12778]';      
 end
if (alpha < alphavec(1)) || (alpha > alphavec(end))
    error('alpha out of range in function kalpha2 !');
end
k0 = interp1(alphavec,kvec,alpha,'linear');
g0 = fcng2(alpha,k0);
k1 = k0*1.001;
if k1 >= kvec(end)
    k1 = k0/1.001;
end
g1 = fcng2(alpha,k1);
itmax = 10; epsi = 1.0e-7; dk = 1.0; it = 0;
while (abs(dk) > epsi) && (it < itmax)
    it = it + 1;
    dk = -g1*(k1 - k0)/(g1 - g0);
    k = k1 + dk;
    k0 = k1;
    k1 = k;
    g0 = g1;
    g1 = fcng2(alpha,k1);
end
function val = fcng2(alfa,k)
% Beregner funksjonen
% g = 1 - sn(alfa,k)
% der sn er en Jacobi-elliptiske funksjon

tol = 1.0e-10; m = k^2;
sn = ellipj(alfa,m,tol);
val = 1 - sn;



