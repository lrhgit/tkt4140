function sol = measles;
solinit = bvpinit(linspace(0,1,5),[0.01; 0.01; 0.01]);
sol = bvp4c(@odes,@bcs,solinit);
x = linspace(0,1,100);
y = deval(sol,x);
plot(x,y(1,:)-0.07,'-',x,y(2,:),'--',x,y(3,:),'-.');
legend('y_1(t) - 0.07  ','y_2(t)','y_3(t)',4);
%=====================================================
function dydx = odes(x,y);
beta = 1575*(1 + cos(2*pi*x));
dydx = [    0.02 - beta*y(1)*y(3)
         beta*y(1)*y(3) - y(2)/0.0279
         y(2)/0.0279 - y(3)/0.01       ];
 function res = bcs(ya,yb);
 res = ya - yb;