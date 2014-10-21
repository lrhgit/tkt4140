% ==== program burger ====
% Beregner fig 4.2 i kapittel 4
%
clear; clf;
x = linspace(-1,1);
x1 = linspace(-1,1.25);
x2 = linspace(1,1.25);
u0 = 1 - x.*x; % t = 0
u1 = 2*(x - 1 + sqrt(2*(1 - x)));  % t = 0.5
u21 = x1 + (sqrt(5 - 4*x1) - 1)/2; % t = 1.0
u22 = x2 - (sqrt(5 - 4*x2) + 1)/2; % t = 1.0
axis([-1 1.5 0 1.2]);
FS = 'FontSize';
hold on
plot(x,u0);
plot(x,u1);
plot(x1,u21);
plot(x2,u22);
hold off
grid
xlabel('x',FS,14);
ylabel('u',FS,14,'Rotation',0);
title('Inviscid Burgers equation',FS,14);
shg;