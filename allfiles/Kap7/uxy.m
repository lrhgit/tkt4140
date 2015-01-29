function uval = uxy(x,y)
% Computes the analytical solution of the Poisson  
% equation  Uxx + Uyy = - 1 on the unit square     
% with u = 0 along the boundaries.

xb = abs(x - 0.5);
yb = abs(y - 0.5);
 if (yb > xb) 
    temp = xb;
    xb = yb;
    yb = temp;
end
zeta = 0.5 - xb;
k = 0;
s = 0.0;
while 1
    k = k + 1;
    n = 2*k -1;
    pin = pi*n;
    t =  cosh(pin*yb)/cosh(pin*0.5)/(n^3);
    ds = t*sin(pin*zeta);
    term = s + t;
    if (k > 200 )|(term == s) 
        uval = (1.0 - zeta)*zeta*0.5 - s*4.0/pi^3;
        return
    end
    s = s + ds;
end
