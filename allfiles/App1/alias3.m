% program alias3
% Plotter punktene xk = k/4, k = 0..8
% Plotter deretter funksjonene f1 = cos(5*Pi*x)
% og f2 = cos(3*Pi*x) som faller sammen i punktene xk.
clear; close
N = 9;
xp = zeros(N,1);
for k = 1:N
    xp(k) = (k - 1)/4;
end
y = cos(5*pi*xp);
x = (0:0.01:2.0)';
f1 = cos(5*pi*x);
f2 = cos(3*pi*x);
axis([-0.2 2.2 -1.2 1.2]);
hold on
plot(xp,y,'ko')
plot(x,f1,'k')
plot(x,f2,'k-.')
hold off
FS = 'FontSize';
title('Aliasing',FS,14)
xlabel('x',FS,14)
ylabel('y',FS,14,'Rotation',0)