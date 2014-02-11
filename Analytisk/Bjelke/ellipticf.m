function value = ellipticf(phi,k)
%
% ellipticf computes the incomplete elliptic integral
% of 1st kind defined by :
%  ellipticf(phi,k) = I1 integrated between 0 and pi/2 where
%       I1 = {[1 + (k*sinu)^2]}^(-1/2)*du
%       Accuracy : Between 13 - 15 digits
%
%       === Other notations ====
%     iellip1(x,kc) = F(phi\90 -alpha)
%              with
%  kc^2 = 1 - k^2 , x = tan(phi)
%  kc = cos(alpha) , k = sin(alpha)
%  m = k^2 , m1 = kc^2 , m + m1 = 1
%
%               === Reference ===
%      W. Press et al. :
%   "Numerical Recipes
%    The Art of Scientific Computing ",
%    Fortran 77, Cambridge University Press 1986
%
tiny = 1.0e-300; ca = 1.49e-8;
x = tan(phi);
kc = sqrt(1 - k^2);
if (abs(x) < tiny)
    value = 0;
    return
end
if (abs(kc) < tiny)
    value = log(tan(pi*0.25 + atan(x)*0.5));
    return
end
qc = kc;
a = 1; b = 1; 
c = x^2; d = 1 + c;
p = sqrt((1 + c*qc^2)/d);
d = x/d; c = 0.5*d/p;
z = a - b;
ey = a;
a = 0.5*(a + b);
y = abs(1/x);
f = 0; l = 0; em = 1;
qc = abs(qc); test = 1;
while test
    b = ey*qc + b; e = em*qc;
    g = e/p; d = f*g + d;
    f = c; ey = a; p = g + p;
    c = 0.5*(d/p + c);
    g = em; em = qc + em;
    a = 0.5*(b/em + a);
    y = -e/y + y;
    if (y == 0)
        y = sqrt(e)*eps;
    end
    test = abs(g - qc) > ca*g;
    if test
        qc = 2.0*sqrt(e);
        l = 2*l;
        if (y < 0)
            l = l + 1;
        end
    end
end
if ( y < 0 )
    l = l + 1;
end
e = (atan(em/y) + pi*l)*a/em;
if (x < 0)
    e = -e;
end
value = e + c*z;      
        