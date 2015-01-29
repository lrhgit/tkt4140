% Program hiemenzplot
% Beregner og plotter fi som funksjon av s= f''(0) ved 
% bruk av skyteteknikk for løsning av Falkner-Skan ligningen
% for beta = 1 (Hiemenz) som funksjon av etainf.
clear
sstart = 0.9; send = 1.3 ; antall = 30;
s = linspace(sstart,send,antall);
fi = s;
options= odeset('RelTol', 1.0e-5);
etainf = 7;
fprintf(' etainf = %6.2f \n',etainf);
xspan = [0 etainf];
for n = 1:antall
   y0 = [0.0 0.0 s(n)];
   [x,y] = ode45('fcnhiemenz',xspan,y0,options);
   fi(n) = y(end,2) - 1.0 ;
end
figure(1)
clf
plot(s,fi)
ymin = min(fi)-0.5;
ylim([ymin 5 ]);
grid
xlabel('s','FontSize',14)
ylabel('\phi','FontSize',14,'FontWeight','Bold')
title('Hiemenz - Nullpunkt for \phi','Fontsize',14,'FontWeight','Bold')