% program four2
% Fourier-utvikling av f(x) = abs(x) for -pi < x < pi
% som en cosinusrekke.
% Plotter først de enkelte modene i fig. 1
% Deretter en 3D-figur av selve funksjonen i fig. 2
clear;close
fac = 4/pi;
x = -pi:0.02:pi;
y = zeros(length(x));
figure(1);
ylim([-1.3 1.8]);
hold on
y = pi/2;
plot(x,y,'k');
for k = 1:2:19
    y = fac*cos(k*x)/k^2;
%     ymax = max(y)
%     ymin = min(y)
    plot(x,y,'k');
end
% axis off
grid
hold off
figure(2);
y = zeros(10,length(x));
z = zeros(size(x));
for k = 1:2:19
    z = z + cos(k*x)/k^2;
    y((k+1)/2,:) = pi/2 - fac*z;
end
plot(y(1:2:9,:)')
% surf(y)
% shading interp
mesh(y)
svart = zeros(64,3);
colormap(svart)
axis off ij
