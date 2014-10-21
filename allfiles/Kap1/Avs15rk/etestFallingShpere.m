% Testprogram for feuler for falling sphere with constant Cd
clear
tspan = [0 1];
y0 = [0 ; 0];
dt = 0.011;
g=9.91;
a=0.007;
[x,y] = feuler('odefun',tspan,y0,dt,g,a);