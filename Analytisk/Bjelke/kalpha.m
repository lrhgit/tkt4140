function k = kalpha(alpha)
% Large deflection of a beam
% Computes k for a given value of alpha
% from the equation k = dn(alpha,k)/cn(alpha,k)/sqrt(2)
% Range : 0 <= alpha <= 5.24636
persistent kvec alphavec
if isempty(kvec)  
    kvec = [0.70711 0.73728 0.76604 0.79335 0.81915 0.84339 0.86603 ...
        0.88701 0.90631 0.92388 0.93969 0.95372 0.96593 0.97630 ...
        0.98481 0.99144 0.99619 0.99657 0.99692 0.99725 0.99756 0.99786 ...
        0.99813 0.99839 0.99863 0.99885 0.99905 0.99923 0.99939 0.99953 ...
        0.99966 0.99976 0.99985 0.99991 0.99996 ]';

    alphavec = [0.00000 0.41836 0.59414 0.73284 0.85475 0.96827 1.07826 ...
          1.18815 1.30088 1.41932 1.54666 1.68689 1.84538 2.03012 ...
          2.25415 2.54156 2.94597 2.99714 3.05107 3.10809 3.16858 ...
          3.23298 3.30183 3.37580 3.45570 3.54258 3.63776 3.74299 ...
          3.86064 3.99406 4.14809 4.33031 4.55336 4.84096 5.24636 ]';
end
if (alpha < 0) | (alpha > alphavec(end))
    error('alpha out of range in function kalpha !');
end
k0 = interp1(alphavec,kvec,alpha,'linear');
g0 = fcng(alpha,k0);
k1 = k0*1.001;
if k1 >= kvec(end)
    k1 = k0/1.001;
end
g1 = fcng(alpha,k1);
itmax = 10; epsi = 1.0e-6; dk = 1.0; it = 0;
while (abs(dk) > epsi) && (it < itmax)
    it = it + 1;
    dk = -g1*(k1 - k0)/(g1 - g0);
    k = k1 + dk;
    k0 = k1;
    k1 = k;
    g0 = g1;
    g1 = fcng(alpha,k1);
end


