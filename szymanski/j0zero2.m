function zj0  = j0zero2(s)
% Function j0zero computes root number s, s = 1,2,.., 
% of the Besselfunction J0 where zj0 is the root,
% using the Newton-Raphsons method.
% Relative error ca. 1.0e-14

% --- Compute an initial value of the root
b0 = (s - 0.25)*pi;
b08 = 0.125/b0;
z0 = b0 + b08*(1 - b08^2*124/3); 
rtol = 1.0e-15;
test = 1.0;
while (test > rtol)
    dz0 = besselj(0,z0)/besselj(1,z0);
    z0 = z0 + dz0;
    test = abs(dz0/z0);
end
zj0 = z0;
