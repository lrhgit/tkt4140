% Program eks2v2eul
% Løser ligning i avsnitt 1.6.2
% med Eulers metode
% Lagrer ikke vektorene
% Ligning : y''(t) = l1*y'(t) + l*y(t) = c,
%           y(0) = 0, y'(0) = 0
%           l1 = -129600, l = 98696, c = 9869.6
clear
l1 = -129600; l = 98696; c = 9869.6;
tmax =  1.2516e-1;
dt = 1.5e-5;
nmax = round(tmax/dt) + 1; % Antall tidskritt
y = 0; v = 0;  t = 0;

% ===  Eulers metode ===
for n = 1 : nmax - 1
    y = y + dt*v ;
    v = (1.0 + l1*dt)*v + dt*(c + l*y);
    t = dt*n;
    fprintf('  %12.4e  %12.5e  %12.5e\n',t,y, v);
end 
fprintf('\n tmax =  %12.4e \n',tmax);
fprintf(' dt =  %12.4e \n',dt);
fprintf(' antall skritt =  %12.0f \n',nmax);
% === Analytisk løsning i endepunktet
a1 = 0.761538735021; a2 = -129600.7615387;
A = 0.09999941239986; B = 5.87600143067797e-7;
ya = A*exp(a1*t) + B*exp(a2*t) - 0.1;
va = a1*A*(exp(a1*t) - exp(a2*t)) ;
fprintf('\n analytisk:  %12.5e  %12.5e \n',ya, va);
erry = abs((y - ya)/ya);
errv = abs((v - va)/va);
fprintf('\n relativ feil:  %10.3e  %10.3e \n',erry, errv);
