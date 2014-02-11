% program kurve
% Finner relasjon mellom theta0 og alpha
% theta0 = t0g i grader
clear
t0g = (0:5:85)';
t0g = [t0g;86;87;88;89];
t0 = t0g*pi/180;
st0 = sin(t0);
k = sqrt(0.5.*(1 + st0));
alf = [0.0 0.41836 0.59414 0.73284 0.85475 0.96827 1.0783 1.1882 1.3009 ...
       1.4193 1.5467 1.6869 1.8454 2.0301 2.2542 2.5416 2.9460 3.6378 ...
       3.8606 4.1481 4.5534 5.2464]';
plot(alf,k)
shg