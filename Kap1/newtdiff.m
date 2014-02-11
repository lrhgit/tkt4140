% Program newtdiff
% Plotting the solution of Newton's
% 1. order equation (1671):
% dy/dx = 1-3*x + y + x^2 +x*y , y(0) = 0
x = linspace(0,2,50)';
a = sqrt(2)/2;
t1 = exp(x.*(1+ x/2));
t2 = erf((1+x)*a)-erf(a);
y = 3*sqrt(2*pi*exp(1))*t1.*t2 + 4*(1-t1)-x;
y1 = x - x.*x + x.^3/3 - x.^4/6;
plot(x,y,'k',x,y1,'k-.')
grid on
FS = 'FontSize';
xlabel('x',FS,14)
ylabel('y','Rotation',0,FS,14)
title('Newton''s equation (1671)',FS,14)