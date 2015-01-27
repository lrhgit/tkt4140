%======== Program Delay34v1 tdma-versjon ==========
% Linearisering etter metoden med etterslep
% Løser eksemplet i avsnitt 3.4
% der diff.ligningen er gitt i ligning 3.4.1
% Bruker Thomas-algoritmen direkte
% Denne versjonen utnytter Matlabs "call by value"-
% teknikk ved kall av TDMA
%==================================================
clear
h = 0.05; % skrittlengde
n = round(1.0/h)-1; % Antall ligninger
fac = (3.0/2.0)*h*h;
nitr = 12; % Antall iterasjoner
%--- Initialisering
%--- a og c blir ikke ødelagt under eliminasjons-prosessen
a = ones(n,1) ; % underdiagonal
c = a ;
ym = zeros(n,1); b = ym; d = ym;
d(n) = - 1.0;
d(1) =  - 4.0;
fprintf('        Itr.      max. avvik  \n');
for k = 1:nitr 
   b = -(2.0 + fac*ym); % hoveddiagonal
   ym1 = tdma(a,b,c,d); % Løser ligningsystemet
   dymax = max(abs((ym1-ym)./ym1));% Beregner relativ avvik
   ym = ym1;
   fprintf(' %10d     %12.3e \n',k,dymax);
 end
 %---- Utskrift av y og feil ----
 xv = [h:h:1.0-h]';
 ya = 4.0./(1 + xv).^2; % Analytisk løsning
 feil = abs((ym1 - ya)./ya);
fprintf('\n      x         y        feil  \n\n')
fprintf('  %7.3f %10.5f  %10.2e \n',[xv ym1 feil]');
