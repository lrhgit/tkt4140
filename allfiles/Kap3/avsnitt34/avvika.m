%========================= program avvika =====================
% Linearisering med bruk av Newton-linearisering
% Lser eksemplet i avsnitt 3.4  der diff.ligningen er gitt
% i ligning 4.1, kap. 3
% Beregner lsning yII som vist i fig. 2.4
% i kap. 2.
% Bruker forskjellig tester for  beregne absolutt avvik
% for et gitt antall iterasjoner.
% Startverdiene beregnes fra en parabel y = -20*x*(1-x)
% som gr gjennom punktene (0,0) og (1,0) og har maksimum 
% y = - 5 for x = 1/2. 
% Bruker tdma for  lse ligningsystemet.
%============================================================
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

h = 0.025; % skrittlengde
ni = 1/h; % Antall intervall
% h m velges slik at ni er et heltall
n = ni-1; % Antall ligninger
a = ones(n,1) ; % underdiagonal
c = a; % overdiagonal 
% a og c blir ikke delagt under eliminasjons-prosessen og 
% kan derfor legges utenfor iterasjonslkka.
x = (h:h:1.0-h)';
ym = -20*x.*(1 - x); % Startverdier
b = zeros(n,1); d = b; dy = b; % allokering
fprintf('        Itr.       \n');

fac = 3.0*h*h;

for it = 1:9 
   b = -(2.0 + fac*ym); % hoveddiagonal
   d = -(fac*0.5)*ym.^2 ; % hyre side 
   d(n) = d(n)- 1.0;
   d(1) = d(1) - 4.0;
   ym1 = tdma(a,b,c,d); % Lser ligningsystemet
   dy = abs(ym1 - ym);
   ta1 = max(dy);
   ta2 = sum(dy)/n;
   ta3 = sqrt(dot(dy,dy))/n;
   ym = ym1; % Oppdatering av y-verdier
   fprintf(' %10d   %9.2e   %9.2e   %9.2e \n',it,ta1, ta2, ta3);
end