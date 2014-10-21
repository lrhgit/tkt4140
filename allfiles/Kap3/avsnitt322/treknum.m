% === treknum===
%
% Lser lign. 2.39 i kompendiet for
% et sammensatt trekantprofil
% iprt > 0 : Utskrift av hvert iprt'te punkt.
% iprt = 0 : Utskrift av frste og siste punkt,
% samt diskontinuitet-verdier.
%

clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

% === Data ===
L = 0.2; D = 0.02; ha = 80; ka = 160; 
hb = 80; kb = 40; 

% === Koeffisienter ===
bta = sqrt(2*ha*L^2/(D*ka));
btb = sqrt(2*hb*L^2/(D*kb));
bt2a = bta*bta; bt2b = btb*btb;

h = 0.01; % skrittlengde
iprt = 5; % printskritt
neq = round(1.0/h);
if mod(neq,2) > eps
    error(' np ikke delelig med 2 ! ');
end
% --- Initialisering ---
a = zeros(neq,1); b = a; c = a; d = a;
theta = a; 

% --- Frste ligning ---
b(1) = 1.0 + 1.5*h;
c(1) = -(1.0 - h*0.5);

m = neq/2 + 1;
for k = 2 : m - 1
    a(k) = -(1.0 -1.5/k);
    b(k) = 2.0*(1.0 - 1.0/k) + bt2a*h/k;
    c(k) = -(1.0 - 0.5/k);
end
for k = m + 1 : neq
    a(k) = -(1.0 - 1.5/k);
    b(k) = 2.0*(1.0 - 1.0/k) + bt2b*h/k;
    c(k) = -(1.0 -0.5/k);
end
a(m) = -4.0*(1.0 - h);
b(m) = 5.0 - 3.0*h + 16.0*h*h;
c(m) = -(1.0 + h);
d(neq) = -c(neq);

theta = tdma(a,b,c,d);
theta = [theta; 1];
% === Beregner theta'(x) fra lign. 2.47) ===
x = [0:h:1]';
n = length(theta);
dtheta = zeros(n,1); % Initialisering
dtheta(1) = bt2a*theta(1);
s = 0;
for k = 2: n
    s = s + 0.5*h*(theta(k) + theta(k-1));
    if k > m
        dtheta(k) = bt2b*s/x(k);
    else
        dtheta(k) = bt2a*s/x(k);
    end
end
fprintf('   x       theta       theta'' \n\n');
x2 = x(m:n);
theta2 = theta(m:n);
dtheta2 = dtheta(m:n);
dtheta2(1) = ka*dtheta(m)/kb;
if iprt > 0
    for k = 1:iprt:m
        fprintf('%6.3f %10.5f  %10.5f \n',x(k),theta(k),dtheta(k));
    end
    if k ~= m
        k = m;
        fprintf('%6.3f %10.5f  %10.5f \n',x(k),theta(k),dtheta(k));
    end
else
    fprintf('%6.3f %10.5f  %10.5f \n',x(1),theta(1),dtheta(1));
    fprintf('%6.3f %10.5f  %10.5f \n',x(m),theta(m),dtheta(m));
end
m2 = n - m + 1;
if iprt > 0
    for k = 1:iprt:m2
        fprintf('%6.3f %10.5f  %10.5f \n',x2(k),theta2(k),dtheta2(k));
    end
    if k ~= m2
        k = m2;
        fprintf('%6.3f %10.5f  %10.5f \n',x2(k),theta2(k),dtheta2(k));
    end
else
    fprintf('%6.3f %10.5f  %10.5f \n',x2(1),theta2(1),dtheta2(1));
    fprintf('%6.3f %10.5f  %10.5f \n',x2(m2),theta(m2),dtheta(m2));
end
    
hp=plot(x,theta);
