% program four1
% Fourier-utvikling av f(x) = 1 for 0 < x < pi
% som en sinusrekke.
% Plotter først de enkelte modene i fig.1
% Deretter en 3d-figur av selve funksjonen.
clear;close;
fac = 4/pi;
x = 0:0.02:pi;
y = zeros(length(x));
figure(1);
xlim([-0.5 3.5]);
ylim([-0.5 1.5]);
hold on
for k = 1:2:19
    y = fac*sin(k*x)/k;
    plot(x,y,'k')
end
% axis off
grid
hold off
figure(2);
y = zeros(10,length(x));
z = zeros(size(x));
for k = 1:2:19
    z = z + sin(k*x)/k;
    y((k+1)/2,:) = fac*z;
end
plot(y(1:2:9,:)')
% surf(y)
% shading interp
mesh(y)
svart = zeros(64,3);
colormap(svart)
axis off ij
