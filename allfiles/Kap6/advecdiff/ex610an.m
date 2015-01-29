%program ex610an
% Beregner den analytiske løsningen av
% "flosshatt" (Top hat) for adveksjon-diffusjons-
% ligningen i avsnitt 6.10, se figur 6.20 - 6.21
% 
clear; close
np = 200; % Antall punkt
np1 = np + 1;
dx = 1/np;
x = (0:dx:1)';
ua = zeros(np1,1); ui = ua;
jc = round(np/2) + 1; % Midtpunkt
nc = jc - 1 ;
hw = round(np*0.1);
jstart = jc - hw;
jend = jc + hw;
for j = jstart: jend
    ui(j) = 0.5;
end
C = 0.5; % Courant-tall
Rc = 5.0; % Reynoldcelletall
D = C/Rc; % Diffusjonstall
n = 0; % Antall tidskritt
antall = 500;
tol = 5.0e-5;
for k = 1 : np
    s = 0;
    j = 0;
    term = 1;
    while (abs(term) > tol) | (j < antall)
        j = j + 1;
        beta = 2*pi*j/np;
        term = exp(-beta^2*D*n)*sin(beta*hw)/j;
        s = s + term*cos(beta*(nc - (k - 1 - C*n)));
    end
    ua(k) = hw/np + s/pi;
end
ua(np1) = ua(1);
plot(x,ua,x,ui);
grid
%     plot(x,u)
%     hold on
%     %axis([-5 5 -1 3])
%     %tell = tell + 1;
%     grid
% 
%         
