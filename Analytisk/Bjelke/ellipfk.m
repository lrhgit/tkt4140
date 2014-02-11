function [K,F] = ellipfk(phi,k)
%
% ellipfk computes the incomplete elliptic integral F
% and the complete integral K of 1st kind defined by :
% 
% 1:  F = I1 integrated between 0 and phi where
%       I1 = {[1 + (k*sinu)^2]}^(-1/2)*du
%              0<= phi <= pi/2
% 2 : K = I1 integrated between 0 and pi/2 where
%
% phi in radians and k^2 < 1.
%
%       Accuracy : app. 13-14 digits
%
%  ===== Reference =====
%    William T. Thompson :
%   "Atlas for Computing Mathematical Functions"
%    John Wiley & Sons, 1997
%
pid2 = pi/2;
a = 1.0; c = 1.0; temp = 1.0;
b = sqrt(1.0 - k^2);
epsd2 = 0.5e-13;
psi = phi;
n = 0;
while abs(c) > epsd2*temp
    tanpsi = tan(psi);
    d = atan((a-b)*tanpsi/(a + b*tanpsi^2));
    psi = 2.0*psi - d;
    temp = a;
    a = 0.5*(a + b);
    c = 0.5*(temp - b);
    b = sqrt(temp*b);
    n = n + 1;
end
K = pid2/a;
F = psi/(2^n*a);
        