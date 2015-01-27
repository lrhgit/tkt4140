% ============== program analyt =========================
% Plotter den analytiske løsningen for en plan dyse.
% f(x) = 2*a*tanh(a*x)
% f'(x) = 2*a^2/(cosh(a*x))^2
% f''(x) = -f(x)*f'(x)
% =======================================================
clear
xinf = 4.0;
a = 2;
x = linspace(0,xinf,50);
ax = a*x;

fx = 2*a*tanh(ax);
fdx = 2*a^2./(cosh(ax).^2);
fd2x = -fx.*fdx;

%st = sprintf('Exercise 5.1 :   p{_1}/p{_0} = % 5.2f', p1p0);
clf
plot(x,fx,'k',x,fdx,'k',x,fd2x,'k')
grid
% FS = 'FontSize'; FW = 'FontWeight';
% xlabel('s',FS,14)
% ylabel('\phi','Rotation',0,FS,14,FW,'Bold')
%title(st,FS,13)
shg