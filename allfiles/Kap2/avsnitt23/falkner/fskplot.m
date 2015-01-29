function fskplot 
% Beregner og plotter fi som funksjon av s = f''(0) ved 
% bruk av skyteteknikk for løsning av Falkner-Skan ligningen
% for en gitt verdi av beta og etainf.
% Bruker nøstede funksjoner
beta = 1;
sstart = 0; send = 1.3 ; antall = 40;
etainf = 4.0;
s = linspace(sstart,send,antall);
fi = s;
options= odeset('RelTol', 1.0e-5);
etaspan = [0 etainf];
for n = 1:antall
   y0 = [0.0 0.0 s(n)];
   [eta,y] = ode45(@fcnfsk,etaspan,y0,options);
   fi(n) = y(end,2) - 1.0 ;
end
clf;
ylow = min(fi)- 0.5;
ylim([ylow 3]);
plot(s,fi)
ylim([ylow 3]);
grid
FS = 'FontSize'; FW = 'FontWeight';
xlabel('s',FS,14)
ylabel('\phi',FS,14,FW,'Bold')
%st = sprintf('\\phi - funksjon. {\\eta}_{\infty} = %5.2f',etainf);
st = sprintf('FSK: \\beta = %4.2f , \\eta_{\\infty} = %5.2f', beta,etainf);
title(st,FS,14,FW,'Bold')
shg;
%--------------------------------------------
function yd = fcnfsk(x,y)
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = y(3);
yd(3) = -(y(1)*y(3) + beta*(1 - y(2)^2));
end
end