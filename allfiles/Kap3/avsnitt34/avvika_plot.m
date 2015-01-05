%========================= program avvika =====================
% Linearisering med bruk av Newton-linearisering
% Loeser eksemplet i avsnitt 3.4  der diff.ligningen er gitt
% i ligning 4.1, kap. 3
%
% Startverdiene beregnes fra en parabel y = -20*x*(1-x)
% som gr gjennom punktene (0,0) og (1,0) og har maksimum 
% y = - 5 for x = 1/2. 
% Initialisering gir loesning yII som vist i fig. 2.4
%
% Bruker forskjellig tester for  beregne absolutt avvik
% for et gitt antall iterasjoner.
% Bruker tdma for  loese ligningsystemet.
%============================================================
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

h = 0.025; % skrittlengde
ni = 1/h; % Antall intervall
% h m velges slik at ni er et heltall

n = ni-1; % Antall ligninger
a = ones(n,1) ; % underdiagonal
c = a; % overdiagonal 

% a og c blir ikke oedelagt under eliminasjons-prosessen og 
% kan derfor legges utenfor iterasjonslkka.

x = (h:h:1.0-h)';
b = zeros(n,1); d = b; dy = b; % allokering

ym = -20*x.*(1 - x);   % Startverdier


fprintf('        Itr.       \n');
fprintf('%10s %11s %11s %11s \n',' ', 'ta1','ta2','ta3');
%fprintf('%10s %11s %11s\n','','tr2','tr3');



% Plot initial solution
 plot(x,ym);
 hold on
 set(gca,'ylim',[-30,10]); % Sett fast y-skala for sammenlikning. 

fac = 3.0*h*h;
 
for it = 1:9 
   b = -(2.0 + fac*ym); % hoveddiagonal
   d = -(fac*0.5)*ym.^2 ; % hyre side 
   d(n) = d(n)- 1.0;
   d(1) = d(1) - 4.0;
   ym1 = tdma(a,b,c,d); % Lser ligningsystemet
   dy = abs(ym1 - ym);
   
   % Absolutte stoppkrit
   ta1 = max(dy);
   ta2 = sum(dy)/n;
   ta3 = sqrt(dot(dy,dy))/n;
   
   %% Relative stoppkrit
   tr2 = sum(dy)/sum(abs(ym1));
   tr3 = max(dy)/max(abs(ym1));
   
   ym = ym1; % Oppdatering av y-verdier
 fprintf(' %10d   %9.2e   %9.2e   %9.2e \n', it, ta1, ta2, ta3);
  %fprintf('%10d   %9.2e   %9.2e\n', it, tr2, tr3);
   pause(1);
   plot(x,ym)
end