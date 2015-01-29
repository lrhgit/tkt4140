% Program eks2v1eul
% Løser eksempel 2 i avsnitt 1.6.1
% med bruk av Eulers metode
% Ligning : y''(t)= l1*y'(t) + l*y(t) + c,
%           y(0) = 0, y'(0) = 0
%           l1 = -129600, l = 98696, c = 9869.6
clear
l1 = -129600; l = 98696; c = 9869.6;
tmax =  1.2516e-2;
dt = 1.5e-5;
nmax = round(tmax/dt) + 1; % Number of time-steps
% Initializing of av vectors
t = zeros(nmax,1); y = t; v = t;
y(1) = 0; v(1) = 0;  t(1) = 0;
% clf;
% === Euler's method ===
for n = 1 : nmax - 1
    y(n+1) = y(n) + dt*v(n);
    v(n+1) = (1.0 + l1*dt)*v(n) + dt*(c + l*y(n));
    t(n+1) = dt*n;
end 
%fprintf('  %12.4e  %12.5e  %12.5e\n',[t y v]');
plot(t,y,'k')
shg
%     hold on
% 
%grid on
%FS = 'FontSize';
%st = sprintf('\\Deltax = %5.2e',dx);
%title(st,FS,13);
%xlabel('x',FS,14);
%ylabel('y',FS,13);