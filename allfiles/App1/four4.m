% program four4
% Fourier-utvikling av f(x) = x for 0 <= x <= pi
% f(x) = -x + 2*pi or pi <= x <= 2*pi
% som en cosinusrekke (Som i tilfelle 2).
clear;close
fac = 4/pi;
x = 0:0.02:2*pi;
% y = pi/2 - fac*cos(x);
% plot(x,y)
% pause
% y = cos(x) + cos(3*x)/3^2;
% y = pi/2 - fac*y;
% plot(x,y)
% pause
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
