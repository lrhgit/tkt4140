% Testprogram for feuler
clear
xs = 1;
xspan = [xs 0];
xd = 2*exp(xs)/sqrt(pi);
ys = erf(xs);
y0 = [ys ; xd];
dx = 0.10;
[x,y] = feuler('odefun2',xspan,y0,dx);
plot(x,y)
grid
shg

