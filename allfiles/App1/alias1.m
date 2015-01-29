% program alias1
% Plotter punktene xk = 0.4*k, k = 0..5
% Plotter deretter funksjonene f1 = cos(3*Pi*x)
% og f2 = cos(2*Pi*x) som faller sammen i punktene xk.
clear; close
N = 6;
xp = zeros(N,1);
for k = 1:6
    xp(k) = (k - 1)*0.4;
end
y = cos(3*pi*xp);
x = (0:0.02:2.0)';
f1 = cos(3*pi*x);
f2 = cos(2*pi*x);
axis([-0.2 2.2 -1.2 1.2]);
hold on
plot(xp,y,'ko')
plot(x,f1,'k')
plot(x,f2,'k-.')
hold off