clear all
clc
close all
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
ax=-MG*x/r^3;
ay=-MG*y/r^3;
i=0;
for t=0:dt:tmax
    i=i+1;
    rp=sqrt(x^2+y^2);
    ax=-MG*x/rp^3;
    ay=-MG*y/rp^3;
    
    up=u+dt*ax;
    vp=v+dt*ay;
    xp=x+dt*u;
    yp=y+dt*v;
    
    x=x+dt*(u+up)/2;
    y=y+dt*(v+vp)/2;
    
    r=sqrt(x^2+y^2);
    u=u+dt*(ax-MG*x/r^3)/2;
    v=v+dt*(ay-MG*y/r^3)/2;
%    u=u_ny;v=v_ny;
    set(h,'XData',x,'YData',y)
    if mod(i,100)==0
        plot(x,y,'b.') 
    end
    drawnow;
end
    