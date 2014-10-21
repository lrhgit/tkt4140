%========================= program fig343 =====================
% Plotter figur i avsnitt 3.4.3
% Linearisering med bruk av Taylorutvikling
% Løser eksemplet i avsnitt 3.4  der diff.ligningen er gitt
% i ligning 4.1, kap. 3
% Beregner løsning yII som vist i fig. 2.4
% i kap. 2.
% Startverdiene beregnes fra en parabel y = -20*x*(1-x)
% som går gjennom punktene (0,0) og (1,0) og har maksimum 
% y = - 5 for x = 1/2. 
% Bruker tdma for å løse ligningsystemet.
%============================================================
clear; clf;
h = 0.025; % skrittlengde
ni = 1/h; % Antall intervall
% h må velges slik at ni er et heltall
n = ni-1; % Antall ligninger
fac = 3.0*h*h;
a = ones(n,1) ; % underdiagonal
c = a; % overdiagonal 
% a og c blir ikke ødelagt under eliminasjons-prosessen og 
% kan derfor legges utenfor iterasjonsløkka.
x = (h:h:1.0-h)';
ym = -20*x.*(1 - x); % Startverdier
b = zeros(n,1); d = b; dy = b; y = ym; % allokering
for it = 1:9 
   b = -(2.0 + fac*ym); % hoveddiagonal
   d = -(fac*0.5)*ym.^2 ; % høyre side 
   d(n) = d(n)- 1.0;
   d(1) = d(1) - 4.0;
   ym1 = tdma(a,b,c,d); % Løser ligningsystemet
   ym = ym1; % Oppdatering av y-verdier
end
x = [0;x;1];
y = [0;y;0];
ym1 = [4;ym1;1];
plot(x,y,'k',x,ym1,'k');
grid
xlabel('x','FontSize',14)
ylabel('y{_s} , y{_{II}}','FontSize',14)
shg;
%---- Utskrift av y ----
% fprintf('\n      x         y   \n\n')
% fprintf('  %7.3f %10.5f \n',[x ym1]');
