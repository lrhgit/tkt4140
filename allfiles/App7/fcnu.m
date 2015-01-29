function uval = fcnu(y,t)
% Computes the analytical solution u(y,t) 
% of the non-steady Couette flow 
% (or heat conduction problem) in appendix 7 
% of the compendium.
% Equation : du/dt = u''(y,t), 0 <= y <= 1
% Initial values : u(y,t) = 0 , t < 0
% Boundary values : u(0,t) = 1 , u(1,t) = 0
ta = abs(t); ya = abs(y);
if (ya > 1|| ya < 0)
    error('y is out of range');
end
if (ta < 0.06)
    % === erfc-serie ===
    arg1 = y/(2*sqrt(t));
    arg2 = (2 - y)/(2*sqrt(t));
    term1 = erfc(arg1);
    term2 = erfc(arg2);
    uval = term1 - term2;
    return
end
% === Sine-serie ===
epsi = 1.0e-10;
term = 1.0;
s = 0.0;
n = 0;
while (term > epsi)
    n = n + 1;
    fac = n*pi;
    term = exp(-fac^2*t)/n;
    s = s + term*sin(fac*y);
end
    uval = 1 - y - 2*s/pi;
    


