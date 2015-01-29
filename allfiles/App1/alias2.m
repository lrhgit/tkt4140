% program alias2
% Plotter punktene xk = k/6, k = 0..12
% Plotter deretter funksjonene f1 = cos(7*Pi*x)
% og f2 = cos(5*Pi*x) som faller sammen i punktene xk.
clear; close
N = 13;
xp = zeros(N,1);
for k = 1:N
    xp(k) = (k - 1)/6;
end
y = cos(7*pi*xp);
x = (0:0.01:2.0)';
f1 = cos(7*pi*x);
f2 = cos(5*pi*x);
axis([-0.2 2.2 -1.2 1.2]);
hold on
plot(xp,y,'ko')
plot(x,f1,'k')
plot(x,f2,'k-.')
hold off