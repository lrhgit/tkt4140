% Testprogram for feuler
clear
xspan = [0 1];
y0 = [4 ; -8];
dx = 0.011;
[x,y] = feuler('odefun',xspan,y0,dx);
