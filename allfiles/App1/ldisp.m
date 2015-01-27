% program ldisp
% Lineær,dispersiv diff. ligning :
%        d(u(x,t)/dt + b^2*d^3(u(x,t))/dt^3 = 0
% Startprofil: u(x,0) = sin(x/b), 0 <= x <= 2*pi
% Løsning : u(x,t) = sin((x + t)/b)
%
% Plotter løsningen for t = 0.3, 0.7 og 1
% Deretter en 3d-figur av selve funksjonen.
clear;close;
x = 0:0.02:2*pi;
y = zeros(length(x));
b = 0.2;
figure(1);
% xlim([-0.5 3.5]);
% ylim([-0.5 1.5]);
hold on
for t = [0.0 5.0]
    y = sin((x + t)/b);
    plot(x,y,'k')
end
%axis off
%grid
hold off
% figure(2);
% y = zeros(10,length(x));
% z = zeros(size(x));
% for k = 1:2:19
%     z = z + sin(k*x)/k;
%     y((k+1)/2,:) = fac*z;
% end
% plot(y(1:2:9,:)')
% % surf(y)
% % shading interp
% mesh(y)
% svart = zeros(64,3);
% colormap(svart)
% axis off ij
