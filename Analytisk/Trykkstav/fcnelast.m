function [x,y] = fcnelast(theta0g,n)
% Trykkstav med stor utbøyning.
% For en gitt verdi av theta0g (i grader)
% beregnes koordinatene (x,y) for den deformerte staven
% Staven deles i n deler slik at (x(1)=0,y(1)=0) er origo i A
% og det siste punktet (x(n+1),y(n+1)) er stavenden B.
% Buelenden l starter i A og ender med l = 1 i B.
t0 = theta0g*pi/180; % i radianer
k = sin(t0/2); k2 = k^2; tol = 1.0e-8;
x = zeros(n+1,1); y = x;
alfa = ellipke(k^2);
for m = 1:n
    l = m/n;
    arg = alfa*l;
    sn = ellipj(arg,k2,tol);
    phi = asin(sn);
    y(m+1) = 2*k*(1-cos(phi))/alfa;
    [Ec,E] = ellipek(phi,k);
    x(m+1) = 2*E/alfa - l;
end
