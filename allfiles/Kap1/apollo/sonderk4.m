% program sonderk4
% Bruker RK4 istedenfor ODE45
clear
global mu lam;
mu = 1/82.45; % R. Newton og Fehlberg
lam = 1 - mu;
%mu = 1/82.2845; % Nyere data fra NASA

% Data fra artikkel av Robert E. Newton
% fprintf('=== Newton-data === \n\n')
% x0 = 1.2; y0 = 0;
% xp0 = 0; yp0 = -1.049;
% tend = 6.192;

% Overstående data modifisert av Fehlberg
fprintf('=== Fehlberg-data === \n\n');
x0 = 1.2; y0 = 0;
xp0 = 0; yp0 = -1.049357509830343;
tend = 6.192169331319632;
dt = 5.6e-4; % tidskritt
steps = fix(tend/dt); % Antall like skritt
ldt = rem(tend,dt); % det siste skrittet.
fprintf(' Tidskritt dt  = %15.5e \n',dt);
fprintf(' Antall skritt steps  = %8.0f \n',steps);
fprintf(' Siste skritt ldt  = %15.8e \n',ldt);

y = [x0 ; y0 ; xp0 ; yp0]; % startverdier
tic
for k = 1:steps
    t = k*dt;
    ynew = RK4C('fcnsonde',t,y,dt);
    y = ynew;
end
% Det siste skrittet
t = t + ldt;
ynew = RK4C('fcnsonde',t,y,ldt);
toc
y = ynew;
xend = y(1); yend =  y(2);
xpend = y(3); ypend =  y(4);

fprintf('\n --- Sluttverdier ---\n\n');
fprintf(' xend  = %15.8e . Avvik: dxend  = %12.4e \n',xend, xend - x0);
fprintf(' yend  = %15.8e . Avvik: dyend  = %12.4e \n',yend, yend - y0);
fprintf(' xpend = %15.8e . Avvik: dxpend = %12.4e \n',xpend, xpend - xp0);
fprintf(' ypend = %15.8e . Avvik: dypend = %12.4e \n',ypend, ypend - yp0);


% Kontroll av Jacobi-integral
fprintf('\n --- Kontroll av Jacobi-integral --- \n\n');
r1 = sqrt((x0 + mu)^2 + y0^2);
r2 = sqrt((x0 - lam)^2 + y0^2);
% Verdi ved t = 0
C = x0^2 + y0^2 + 2*lam/r1 + 2*mu/r2 - (xp0^2 + yp0^2);
fprintf('Ved start : C = %12.6e \n',C);
% Verdi ved tend
r1 = sqrt((xend + mu)^2 + yend^2);
r2 = sqrt((xend - lam)^2 + yend^2);
Cs = xend^2 + yend^2 + 2*lam/r1 + 2*mu/r2 - (xpend^2 + ypend^2);
fprintf('Ved slutt : C = %12.6e \n',Cs);
fprintf('Avvik: dC = %12.6e \n',Cs - C);