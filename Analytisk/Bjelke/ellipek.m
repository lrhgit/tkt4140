function [Ec,E] = ellipek(phi,k)
%
% ellipek computes the incomplete elliptic integral E
% and the complete integral Ec of 2st kind defined by :
% 
% 1:  E = I1 integrated between 0 and phi where
%       I1 = sqrt(1 + (k*sinu)^2)*du
%              0<= phi <= pi/2
% 2 : Ec = I1 integrated between 0 and pi/2.
%
%      phi is in radians and k^2 <= 1.
%
%       Accuracy : app. 13 - 14 digits
%
%  ===== Reference =====
%    William T. Thompson :
%   "Atlas for Computing Mathematical Functions"
%    John Wiley & Sons, 1997
%
pid2 = pi/2;
m = k^2;
if ( m == 1.0)
    Ec = 1.0;
    E = sin(phi);
    return
end
a = 1.0; c = 1.0; temp = 1.0;
b = sqrt(1.0 - m);
eps2 = 0.5e-13;
psi = phi;
n = 0;
sum1 = 1.0 - 0.5*m;
sum2 = 0.0;
power2 = 1.0;
while abs(c) > eps2*temp
    tanpsi = tan(psi);
    d = atan((a-b)*tanpsi/(a + b*tanpsi^2));
    psi = 2.0*psi - d;
    c = 0.5*(a - b);
    temp = a;
    a = 0.5*(a + b);
    c = 0.5*(temp - b);
    b = sqrt(temp*b);
    sum1 = sum1 - power2*c^2;
    sum2 = sum2 + c*sin(psi);
    power2 = power2*2.0;
    n = n + 1;
end
Ec = sum1*pid2/a;
E = sum2 + sum1*psi/(2^n*a);
        