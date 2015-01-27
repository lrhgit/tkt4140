% Testprogram for heun
clear
xspan = [0 1];
y0 = [4 ; -8];
dx = 0.1;
[x,y] = heun('odefun',xspan,y0,dx);
