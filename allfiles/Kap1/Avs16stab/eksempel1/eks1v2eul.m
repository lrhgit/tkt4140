% Program eks1v2eul
% Løser eksempel 1 i avsnitt 1.6.1
% med Eulers metode
% Lagrer ikke vektorene

% Ligning : x''(t) = l1*x'(t) + l*x(t) 
%           x(0) = 1, x'(0) = 0 = v(0)
%           l1 = -2*c, l = -1 
%           Bruker  c = 100 -> l1 =-200

clear
c = 100;
l1 = -2*c; l = -1; 
tmax =  20;
dt = 0.009;
nmax = round(tmax/dt) + 1; % Antall tidskritt
x = 1; v = 0;  t = 0;

% ===  Eulers metode ===
for n = 1 : nmax - 1
    x = x + dt*v ;
    v = (1.0 + l1*dt)*v + dt*l*x;
    t = dt*n;
    fprintf('  %12.4e  %12.5e  %12.5e\n',t,x, v);
end 
fprintf('\n tmax =  %12.4e \n',tmax);
fprintf(' dt =  %12.4e \n',dt);
fprintf(' antall skritt =  %12.0f \n',nmax);
% === Analytisk løsning i endepunktet
d = sqrt(c^2 - 1);
b1 = -c + d; b2 = -c -d; b12 = 2*d;
A = -b2/b12; B = b1/b12;
xa = A*exp(b1*t) + B*exp(b2*t) ;
va = b1*A*exp(b1*t) + b2*B*exp(b2*t) ;
fprintf('\n analytisk:  %12.5e  %12.5e \n',xa, xa);
errx = abs((x - xa)/xa);
errv = abs((v - va)/va);
fprintf('\n relativ feil:  %10.3e  %10.3e \n',errx, errv);
