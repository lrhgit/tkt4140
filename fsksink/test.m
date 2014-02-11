%program test
a = 2;
x = linspace(0,2*a);
f = sign(x-a).*sqrt(abs(x-a));
plot(x,f)
grid on