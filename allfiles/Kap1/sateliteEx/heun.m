clear all; close all; clc;
FS = 20; set(0,'DefaultLineLineWidth',3,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);
global MG
MG=1;
x=1;
y=0;
u=0;
v=0.9;
figure(1)
plot(0,0,'ro','MarkerSize',40)
Xmin=-1.5; Xmax=-Xmin; Ymin=Xmin; Ymax=-Xmin; axis([Xmin Xmax Ymin Ymax])
hold on
h=plot(x,y,'bo','MarkerSize',20);
tmax=25;
dt=0.01;
i=0;

r=sqrt(x^2+y^2);

E0=.5*(u^2+v^2)-MG/r;
H0=x*v-y*u;

for t=0:dt:tmax
    i=i+1;
    
    xn = x;
    yn = y;    

    up=u+dt*a(xn,yn);
    vp=v+dt*a(yn,xn);

    xp=x+dt*up;
    yp=y+dt*vp;

    x=x+dt*(u+up)/2;
    y=y+dt*(v+vp)/2;
    
    u=u+dt*(a(xn,yn)+a(xp,yp))/2;
    v=v+dt*(a(yn,xn)+a(yp,xp))/2;

    set(h,'XData',x,'YData',y)
     if mod(i,100)==0
         plot(x,y,'b.')
          r=sqrt(x^2+y^2);
        E=.5*(u^2+v^2)-MG/r;    
        display(sprintf('Energy error: %6.4e',(E0-E)/E0));
        H=x*v-y*u;
        display(sprintf('Spinn error: %6.4e',(H0-H)/H0));
     end
    drawnow;
    
end
    
