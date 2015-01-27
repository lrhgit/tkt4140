% Program beuler2
% Løser eksempel 2 i avsnitt 1.6.1
% med baklengs Eulers metode
% Lagrer ikke vektorene
% Ligning : y''(t)= l1*y'(t) + l*y(t) + c,
%           y(0) = 0, y'(0) = 0
%           l1 = -129600, l = 98696, c = 9869.6
clear
l1 = -129600; l = 98696; c = 9869.6;
tmax =  0.12516;
dt = 6.0e-5;
nmax = round(tmax/dt) + 1; % Antall tidskritt
y = 0; v = 0;  t = 0;
fac = 1/(1 + dt*(-l1 - l*dt));
% === Baklengs Eulers metode ===
for n = 1 : nmax - 1
    y = fac*((1 - dt*l1)*y + dt*(v + dt*c));
    v = fac*(v + dt*l*y + dt*c);
    t = dt*n;
    fprintf('  %12.4e  %12.5e  %12.5e\n',t,y, v);
end 
fprintf('\n tmax =  %12.4e \n',tmax);
fprintf(' dt =  %12.4e \n',dt);
% === Analytisk løsning i endepunktet
a1 = 0.761538735021; a2 = -129600.7615387;
A = 0.09999941239986; B = 5.87600143067797e-7;
ya = A*exp(a1*t) + B*exp(a2*t) - 0.1;
va = a1*A*exp(a1*t) + a2*B*exp(a2*t) ;
fprintf('\n analytisk:  %12.5e  %12.5e \n',ya, va);
%fprintf('  %12.4e  %12.5e  %12.5e\n',[t y v]');
%      plot(t,y,'k')
%     hold on
% for dx = [ 0.3 1.0 ]
%     nmax = round(xmax/dx); % Number of time-steps
%     % Initializing of av vectors
%     x = zeros(nmax,1); y = x; 
%     y(1) = 1;  x(1) = 0;
%     % === Euler's method ===
%     for n = 1 : nmax - 1
%         y(n+1) = (1 -  a*dx)*y(n);
%         x(n+1) = dx*n;
%     end 
%     plot(x,y,'k')
% end
%grid on
%FS = 'FontSize';
%st = sprintf('\\Deltax = %5.2e',dx);
%title(st,FS,13);
%xlabel('x',FS,14);
%ylabel('y',FS,13);