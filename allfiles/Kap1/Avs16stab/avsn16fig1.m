x = linspace(0,1);
l1= 1;
y1 = exp(l1*x);
l2 = 2.0;
y2 = exp(-l2*x);
plot(x,y1,'k',x,y2,'k','LineWidth',1);
xlabel('x','FontSize',14)
ylabel('y','FontSize',14,'Rotation',0)
title('Funksjonen y = exp(\lambdax)','FontSize',14)
grid on;
shg;