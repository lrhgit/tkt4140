function zj0  = j0zero(s)
% Function j0zero computes root number s, s = 1,2,.., 
% of the Besselfunction J0 where zj0 is the root,
% using table-values, asymptotic formulae and 
% Newton-Raphsons method.
% Relative error app. 1.0e-14 to 1.0e-15
%
% --- First 5 zeros from table
%
 sz(1) = 2.404825557695773;
 sz(2) = 5.520078110286311;
 sz(3) = 8.653727912911012;
 sz(4) = 11.79153443901428;
 sz(5) = 14.93091770848779;
if ( s <= 5)
    zj0 = sz(s);
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
% --- One iteration with Newton-Raphson's method
dz0 = besselj(0,z0)/besselj(1,z0);
z0 = z0 + dz0;
zj0 = z0;