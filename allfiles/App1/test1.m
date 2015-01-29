%program test1
clear; close
x = linspace(-4,4);
y = ones(size(x));
k = (x < 0);
y(k) = sin(x(k));
k = (0 < x & x <=1);
y(k) = x(k);
plot(x,y);
ylim([-1.5 1.5]);
grid