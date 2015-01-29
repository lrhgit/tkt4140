clear all; close all; clc;
FS = 20;
set(0,'DefaultLineLineWidth',3,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);
% set(0,'DefaultTextFontSize',FS);
% set(0,'DefaultTextFontName','Helvetica');
MG=1;
x=1;
y=0;
u=0;
v=0.7;
figure(1)
plot(0,0,'ro','MarkerSize',40)
axis([-1.5 1.5 -1.5 1.5])
hold on
h=plot(x,y,'bo','MarkerSize',20);
tmax=10;
dt=0.001;
i=0;
for t=0:dt:tmax
    i=i+1;
    r=sqrt(x^2+y^2);
    ax=-MG*x/r^3;
    ay=-MG*y/r^3;
    u_ny=u+dt*ax;
    v_ny=v+dt*ay;
    x=x+dt*u;
    y=y+dt*v;
    u=u_ny;v=v_ny;
    set(h,'XData',x,'YData',y)
    if mod(i,100)==0
        plot(x,y,'b.') 
    end
    drawnow;
end
    