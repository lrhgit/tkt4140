function zj0  = j0test(s)
% Function j0zero computes root number s, s = 1,2,.., 
% of the Besselfunction J0 where zj0 is the root,
% using table-values, asymptotic formulae and 
% Newton-Raphsons method.
% Relative error app. 1.0e-15

if ( s == 1)
    zj0 = 2.404825557695773;
    return
end

% --- Compute an initial value z0 of the root
% --- from asymptotic formulae 9.5.12
% --- in Handbook of Mathematical Functions

b0 = (s - 0.25)*pi;
b08 = 0.125/b0; b082 = b08^2; 
z0 = b0 + b08*(1 - b082*(124/3 - b082*120928/15)); 

if (s > 30)
    zj0 = z0;
    return
end
% --- Newton-Raphson's method

dz0 = besselj(0,z0)/bessel(1,z0);
z0 = z0 + dz0;
if (s > 5)
    zj0 = z0;
    return
end
dz0 = besselj(0,z0)/bessel(1,z0);
z0 = z0 + dz0;
zj0 = z0;
