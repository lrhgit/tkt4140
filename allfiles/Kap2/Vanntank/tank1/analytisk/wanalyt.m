% program wanalyt
% Avsnitt 2.4.1: Vanntank
% Beregner analytisk løsning av en sirkulær 
% vanntank der veggtykkelsen er konstant.
%
% Beregner og plotter utbøyningen w(x) 
% som funksjon av beta
clear
clf
% R = 9.0; % Tankradius
% H = 8.0; % Høyde
% t = 0.35; % Veggtykkelse 
% ny = 0.2;  % Poissons tall
% 
% b = (3*(1 - ny^2)*H^4/(R*t)^2)^0.25;

bstart = 1; bslutt = 6; bstep = 1;
dx = 0.01; kmax = 1 + round(1/dx);
xv = (0 : dx :1.0)';
w = zeros(kmax,1);
FW = 'FontWeight';
hold on
for b = bstart:bstep:bslutt
    %fprintf('beta = %12.5e \n',b);
    s2b =sin(2*b); cb = cos(b); cb2 = cb^2;
    J = 4*b*(cb2 + cosh(b)^2);
    C1 = (2*cb2*(b-1) - s2b*b + b*(1+exp(-2*b)))/J;
    C2 = (2*cb2*b + s2b*(b-1) - (1+b)*(1+exp(-2*b)))/J;
    C3 = (2*cb2*(1+ b) + s2b*b + b*(1+exp(2*b)))/J;
    C4 = -(2*cb2*b -s2b*(1+ b) + (1-b)*(1+exp(2*b)))/J;
    %     wp og dwp er bidraget fra partikulærløsningen
    for k = 1:kmax 
        x = xv(k);
        wp = -(1 - x);
        sbx = sin(b*x); cbx = cos(b*x);
        wh = exp(b*x)*(C1*cbx + C2*sbx) + exp(-b*x)*(C3*cbx + C4*sbx);
        w(k) = wh + wp;
    end
    plot(xv,w)
end
hold off
title('Utbøyning w(x) som funksjon av {\beta}',FW,'Bold');
xlabel('x',FW,'Bold')
ylabel('w(x)',FW,'Bold')
shg
