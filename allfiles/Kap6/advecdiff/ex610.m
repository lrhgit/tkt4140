%program ex610
% Beregner den numeriske og analytiske løsningen av
% "flosshatt" (Top hat) for adveksjon-diffusjons-
% ligningen i avsnitt 6.10, se figur 6.20 - 6.21
%  --- 0 <= theta <= 1 ---
%  theta = 0 gir eksplisitt metode
%  theta = 1 gir total implisitt metode
% --- 0 <= alpha <= 1 ---
% alpha = 0 gir sentraldifferanser
% alpha = 1 gir full oppstrøms-differanser
clear; close
np = 100; % Antall punkt
np1 = np + 1;
dx = 1/np;
x = (0:dx:1)';
% Allokering
ua = zeros(np1,1); u = ua; ui = ua; 
d = zeros(np1,1); e = d; h = d;
%
jc = round(np/2) + 1; % Midtpunkt
nc = jc - 1 ;
hw = round(np*0.1);
jstart = jc - hw;
jend = jc + hw;
% Startverdier 
for j = jstart:jend
    ui(j) = 0.5;
end
C = 0.5; % Courant-tall
Rc = 2.0; % Reynoldcelletall
D = C/Rc; % Diffusjonstall
%
n = 200; % Antall tidskritt
%
fprintf('Reynoldscelletall Rc = %6.2f \n',Rc);
fprintf('Courant-tall C = %6.2f \n',C);
fprintf('Reynoldscelletall Rc = %6.2f \n',Rc);
fprintf('Antall tidskritt = %5.0f \n',n);
%
% --- Numerisk løsning ---
%
alpha = 0.5; theta = 0;
fprintf('alpha = %6.2f \n',alpha);
fprintf('theta = %6.2f \n',theta);
Rcinv = 1/Rc;
tp1 = C*((1 - alpha)*0.5 - Rcinv);
c = theta*tp1; % overdiagonal
cn = -(1 - theta)*tp1;
tp2 = C*(alpha + 2*Rcinv);
b = 1 + theta*tp2; % hovediagonal
bn = 1 - (1 - theta)*tp2;
tp3 = C*((1 + alpha)*0.5 + Rcinv);
a = -theta*tp3; % underdiagonal
an = (1 - theta)*tp3;
u = ui; % Startverdier
% --- Start løkke ---
cycle = 0;
while cycle <=n
    % Beregn høyre side
    for k = 2:np
        d(k) = cn*u(k+1) + bn*u(k) + an*u(k-1);
    end
    % Startverdier for tdma
    d(1) = (cn*u(2) + bn*u(1) + an*u(np))/b;
    e(1) = - c/b;
    h(1) = -a/b;
    p2 = c;
    p3 = -d(np);
    % Beregn syklus
    for k = 2:np-1
        q = -1/(b + a*e(k - 1));    
        e(k) = c*q;  
        h(k) = a*h(k-1)*q;   
        d(k) = (-d(k) + a*d(k-1))*q;    
        p3 = p3 + p2*d(k-1);    
        p2 = p2*e(k-1);
    end   
    p1 = 1/q;
    temp = (p2 + a)*d(np -1) + p3;
    u(np) = temp/(p1 - (p2 + a)*(e(np - 1) + h(np -1)));
    umax = abs(u(np));
    for j = np -1:-1:1   
        u(j) = e(j)*u(j+1)+h(j)*u(np) + d(j);    
        auj = abs(u(j));    
        if auj > umax        
            umax = auj;    
        end        
    end    
    u(np1) = u(1);
    cycle = cycle + 1;
end
%
% --- Analytisk løsning ---
%
antall = 500;
tol = 5.0e-6;
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
%  --- Plotting ---
% 1) For Rc = 2.0 og alpha = 0:
%    daspect([ 1 1.5 1])
%    ylim([ -0.1 0.5])
% 2) For Rc = 2.0 og alpha = 0.5 og 1.0:
%    daspect([ 1 1 1])
% For Rc = 3.95 og Rc = 5.0:
%    daspect([ 1 3 1])
%
ua(np1) = ua(1);
plot(x,u,x,ua,'k-.');

%daspect([1 1 1]);
grid
%ylim([-0.1 0.5]);
FS = 'FontSize'; FW = 'FontWeight';
% Overskrift for fig. 6.20b
% title('\alpha = 0, \theta = 0, R_c = 2.0',FS,14,FW,'Bold')
% Overskrift for fig. 6.20c
% title('\alpha = 0, \theta = 0, R_c = 3.95',FS,14,FW,'Bold')
% % Overskrift for fig. 6.20d
% title('\alpha = 0, \theta = 0, R_c = 5.0',FS,14,FW,'Bold')
% % Overskrift for fig. 6.21a
% title('\alpha = 1.0, \theta = 0, R_c = 2.0',FS,14,FW,'Bold')
% % Overskrift for fig. 6.21b
title('\alpha = 0.5, \theta = 0, R_c = 2.0',FS,14,FW,'Bold')

