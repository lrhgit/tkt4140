% program sonde
clear global mu lam;
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

time = [0 tend];
ystart = [x0 ; y0 ; xp0 ; yp0];
relfeil = 1.0e-6; absfeil = 1.0e-4*relfeil;
options = odeset('RelTol',relfeil,'AbsTol',absfeil,'Refine',1);
%tic
[t,y] = ode45(@fcnsonde,time,ystart,options);
%toc
% plot(t,y(:,1),t,y(:,2));
% grid
fprintf('Relativ feil = %8.3e \n',relfeil);
fprintf('Antall skritt = %8.0f \n',length(t));
% max-verdier for dt
[maxdt,ind1] = max(diff(t));
fprintf('\n Maksimalverdi for dt: \n');
fprintf(' dt = %8.3e for t = %8.3e \n',maxdt,t(ind1));
fprintf(' x = %8.3e ,  y = %8.3e \n', y(ind1,1),y(ind1,2));
% min-verdier når sonden er nærmest jorda
fprintf('\n Minimalverdi for dt: \n');
r1 = sqrt((y(:,1) + mu).^2 + y(:,2).^2);
[minr, ind2] = min(r1);
mindt = t(ind2+1) - t(ind2);
fprintf(' dt = %8.3e for t = %8.3e \n',mindt,t(ind2));
fprintf(' x = %8.3e ,  y = %8.3e \n', y(ind2,1),y(ind2,2));

xend = y(end,1); yend =  y(end,2);
xpend = y(end,3); ypend =  y(end,4);

fprintf('\n --- Sluttverdier ---\n\n');
fprintf(' xend  = %15.8e . Avvik: dxend  = %12.3e \n',xend, xend - x0);
fprintf(' yend  = %15.8e . Avvik: dyend  = %12.3e \n',yend, yend - y0);
fprintf(' xpend = %15.8e . Avvik: dxpend = %12.3e \n',xpend, xpend - xp0);
fprintf(' ypend = %15.8e . Avvik: dypend = %12.3e \n',ypend, ypend - yp0);

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