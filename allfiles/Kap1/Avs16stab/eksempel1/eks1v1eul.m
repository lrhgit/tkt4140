% Program eks1v1eul
% Løser eksempel 1 i avsnitt 1.6.1
% av Eulers metode.
% Plotter løsning

% Ligning : x''(t)= l1*x'(t) + l*x(t)
%           x(0) = 1, x'(0) = 0 = v(0)
%           l1 = -2*c, l = -1
%           Bruker c = 100 -> l1 = -200
clear
c = 100;
l1 = -2*c; l = -1;
tmax =  0.22;
dt = 0.011;
nmax = round(tmax/dt) + 1; % Number of time-steps
% Initializing of av vectors
t = zeros(nmax,1); x = t; v = t;
x(1) = 1; v(1) = 0;  t(1) = 0;
% clf;
% === Euler's method ===
for n = 1 : nmax - 1
    x(n+1) = x(n) + dt*v(n);
    v(n+1) = (1.0 + l1*dt)*v(n) + dt*l*x(n);
    t(n+1) = dt*n;
end 
%fprintf('  %12.4e  %12.5e  %12.5e\n',[t y v]');
plot(t,x,'k')
shg
%     hold on
% 
%grid on
%FS = 'FontSize';
%st = sprintf('\\Deltax = %5.2e',dx);
%title(st,FS,13);
%xlabel('x',FS,14);
%ylabel('y',FS,13);