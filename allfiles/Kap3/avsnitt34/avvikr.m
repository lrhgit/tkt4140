<<<<<<< HEAD
%========================= program avvikr =====================
% Linearisering med bruk av Taylorutvikling
% Løser eksemplet i avsnitt 3.4  der diff.ligningen er gitt
% i ligning 4.1, kap. 3
% Beregner løsning yII som vist i fig. 2.4
% i kap. 2.
% Bruker forskjellig tester for � beregne relativt avvik
% for et gitt antall iterasjoner.
% Startverdiene beregnes fra en parabel y = -20*x*(1-x)
% som går gjennom punktene (0,0) og (1,0) og har maksimum 
% y = - 5 for x = 1/2. 
% Bruker tdma for å løse ligningsystemet.
%============================================================
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

h = 0.022; % ønsket skrittlengde
ni = round(1/h); % Antall intervall som må være heltall
h = 1/ni;  %Justerer h slik at ni blir heltall
fprintf('Skrittlengde er h: %f\n',h)
n = ni-1; % Antall ligninger

fac = 3.0*h*h;
a = ones(n,1) ; % underdiagonal
c = a; % overdiagonal 
% a og c blir ikke overskrevet i eliminasjons-prosessen og 
% kan derfor legges utenfor iterasjonsløkka.

x = (h:h:1.0-h)';

%ym= yI(x)+100*rand(n,1);
%ym=100*ones(n,1);

%ym = -20*x.*(1 - x); % Startverdier
ym=-100*ones(n,1);

b = zeros(n,1); d = b; dy = b; % allokering
fprintf('        Itr.       \n');
MaxIt = 15; 

for it = 1:MaxIt 
   b = -(2.0 + fac*ym); % hoveddiagonal
   d = -(fac*0.5)*ym.^2 ; % høyre side
   d(1) = d(1) - 4.0;
   d(n) = d(n)- 1.0;
   
   ym1 = tdma(a,b,c,d); % Løser ligningsystemet
   
   dy = abs(ym1 - ym);
   tr3 = max(dy)/max(abs(ym1));
   tr2 = sum(dy)/sum(abs(ym1));
   ym = ym1; % Oppdatering av y-verdier
   fprintf(' %10d   %9.2e   %9.2e  \n',it,tr2, tr3);
end

plot(x,yI(x),x,ym,'--')
xlabel('x')
ylabel('y')
legend('analytical','numerical')

=======
%========================= program avvikr =====================
% Linearisering med bruk av Taylorutvikling
% L??ser eksemplet i avsnitt 3.4  der diff.ligningen er gitt
% i ligning 4.1, kap. 3
% Beregner l??sning yII som vist i fig. 2.4
% i kap. 2.
% Bruker forskjellig tester for ??? beregne relativt avvik
% for et gitt antall iterasjoner.
% Startverdiene beregnes fra en parabel y = -20*x*(1-x)
% som g??r gjennom punktene (0,0) og (1,0) og har maksimum 
% y = - 5 for x = 1/2. 
% Bruker tdma for ?? l??se ligningsystemet.
%============================================================
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

h = 0.022; % ??nsket skrittlengde
ni = round(1/h); % Antall intervall som m?? v??re heltall
h = 1/ni;  %Justerer h slik at ni blir heltall
fprintf('Skrittlengde er h: %f\n',h)
n = ni-1; % Antall ligninger

fac = 3.0*h*h;
a = ones(n,1) ; % underdiagonal
c = a; % overdiagonal 
% a og c blir ikke overskrevet i eliminasjons-prosessen og 
% kan derfor legges utenfor iterasjonsl??kka.

x = (h:h:1.0-h)';

%ym= yI(x)+100*rand(n,1);
%ym=100*ones(n,1);

ym = -20*x.*(1 - x); % Startverdier
%ym=-100*ones(n,1);

b = zeros(n,1); d = b; dy = b; % allokering
fprintf('        Itr.       \n');
MaxIt = 15; 

for it = 1:MaxIt 
   b = -(2.0 + fac*ym); % hoveddiagonal
   d = -(fac*0.5)*ym.^2 ; % h??yre side
   d(1) = d(1) - 4.0;
   d(n) = d(n)- 1.0;
   
   ym1 = tdma(a,b,c,d); % L??ser ligningsystemet
   
   dy = abs(ym1 - ym);
   tr3 = max(dy)/max(abs(ym1));
   tr2 = sum(dy)/sum(abs(ym1));
   ym = ym1; % Oppdatering av y-verdier
   fprintf(' %10d   %9.2e   %9.2e  \n',it,tr2, tr3);
end

plot(x,ym,'--')
xlabel('x')
ylabel('y')
%legend('analytical','numerical')

>>>>>>> master
