%========================= Program delay34v2 =====================
% Linearisering med bruk av metoden med etterslep.
% Løser eksemplet i avsnitt 3.4
% der diff.ligningen er gitt i ligning 4.1, kap. 3
% Forsøker å beregne både løsning yI og yII som vist i fig. 2.4
% i kap. 2. Hensikten er å finne hvor nært vi må tippe 
% startverdiene forat vi skal få konvergens mot de respektive 
% løsningene.
% Startverdiene beregnes fra en parabel y = -4*p*x*(1-x)
% som går gjennom punktene (0,0) og (1,0). 
% y-avstanden fra x-aksen til minimumspunktet betegnes p.
% For p = 0 blir y = 0.
% Ved å velge p = 1,2, osv., vil parabelen nærme seg den 
% virkelige løsningen for yII.
% Bruker tdma for å løse ligningsystemet.
%==================================================================
clear
h = 0.05; % skrittlengde
ni = 1/h; % Antall intervall
% h må velges slik at ni er et heltall
n = ni-1; % Antall ligninger
fac = 1.5*h*h;
a = ones(n,1) ; % underdiagonal
c = a; % overdiagonal 
% a og c blir ikke ødelagt under eliminasjons-prosessen og 
% kan derfor legges utenfor iterasjonsløkka.
p = 20;
x = (h:h:1.0-h)';
ym = -4*p*x.*(1 - x); % Startverdier
b = zeros(n,1); d = b; % allokering
d = 0 ; % høyre side 
d(n) = - 1.0;
d(1) = - 4.0;
it = 0; itmax = 15; dymax = 1.0; RelTol = 1.0e-5;
fprintf('        p = %6.2f \n\n',p);
fprintf('        Itr.      max. avvik  \n');
while (dymax > RelTol) & (it < itmax)
   it = it + 1;
   b = -(2.0 + fac*ym); % hoveddiagonal
   ym1 = tdma(a,b,c,d); % Løser ligningsystemet
   dymax = max(abs((ym1-ym)./ym1));% Beregner relativ avvik
   ym = ym1; % Oppdatering av y-verdier
   fprintf(' %10d     %12.3e \n',it,dymax);
end
%---- Utskrift av y ----
x = [0;x;1];
ym1 = [4;ym1;1];
fprintf('\n      x         y   \n\n')
fprintf('  %7.3f %10.5f \n',[x ym1]');
