function fskplot2 
% Beregner og plotter fi som funksjon av s= f''(0) ved 
% bruk av skyteteknikk for løsning av Falkner-Skan ligningen
% for en gitt verdi av beta og flere verdier av etainf.
clear; clear global beta;
global beta;
beta = 1.0;
sstart = 1.1; send = 1.25 ; antall = 50;
s = linspace(sstart,send,antall);
fi = s;
options= odeset('RelTol', 1.0e-5);
clf;
hold on
for etainf = [8 12]
    etaspan = [0 etainf];
    for n = 1:antall
        y0 = [0.0 0.0 s(n)];
       [eta,y] = ode45(@fcnfsk,etaspan,y0,options);
       fi(n) = y(end,2) - 1.0 ;
    end
    plot(s,fi,'k')
    ylim([-3 3]);
end
hold off
grid
FS = 'FontSize'; FW = 'FontWeight';
xlabel('s',FS,14)
ylabel('\phi',FS,14,FW,'Bold')
st = sprintf('FSK: \\beta = %4.2f ', beta);
title(st,FS,14,FW,'Bold')
shg;
%-------------------------------------------
function yd = fcnfsk(x,y)
% Used by fskplot2
global beta;
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = y(3);
yd(3) = -(y(1)*y(3) + beta*(1 - y(2)^2));