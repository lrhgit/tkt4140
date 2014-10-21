clear all
clc
close all
global MG
MG=1;
x=1;
y=0;
u=0;
v=.7;
figure(1)
plot(0,0,'ro','MarkerSize',40)
axis([-1.5 1.5 -1.5 1.5])
hold on
h=plot(x,y,'bo','MarkerSize',20);
tmax=100;
dt=0.01;
r=sqrt(x^2+y^2);
%ax=-MG*x/r^3;
%ay=-MG*y/r^3;
i=0;
for t=0:dt:tmax
    i=i+1;
    
    u=u+dt*a(x,y);
    v=v+dt*a(y,x);
    
    x=x+dt*u;
    y=y+dt*v;
    
%    u=u_ny;v=v_ny;
    set(h,'XData',x,'YData',y)
    if mod(i,100)==0
        plot(x,y,'b.') 
    end
    drawnow;
end
    
