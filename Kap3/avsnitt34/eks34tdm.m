%======== Program eks34 tdma-versjon =======
% Løser eksemplet i avsnitt 3.4
% der diff.ligningen er gitt i ligning 3.4.1
% Bruker tdma som løser.
%============================================
clear
h = 0.1; % skrittlengde
n = round(1.0/h)-1; % Antall ligninger
fac = 3.0*h*h;
nitr = 6; % Antall iterasjoner
%--- a  blir ikke ødelagt under eliminasjons-prosessen
a = ones(n,1) ; % underdiagonal
%--- Initialisering
c = zeros(n,1); % overdiagonal 
y = c; b = y; d = y; 
fprintf('        Itr.      max. resid  \n');
for k = 1:nitr 
   c = ones(n,1);
   b = -(2.0 + fac*y); % hoveddiagonal
   d = -(fac*0.5)*y.^2 ; % høyre side 
   d(n) = d(n)- 1.0;
   d(1) = d(1) - 4.0;
   %-------------------------------------- 
   % Løser ligningsystemet
   %--------------------------------------
   x = tdma(a,b,c,d); % Løs ligningsystemet
   %--- Beregn relativ avvik
   dymax = max(abs((x-y)./x));
   y = x;
   fprintf(' %10d     %12.3e \n',k,dymax);
 end
 %---- Utskrift av y og dy ----
 xv = [h:h:1.0-h]';
 ex = 4.0./(1 + xv).^2;
 feil = abs((y - ex)./ex);
fprintf('\n      x         y        feil  \n\n')
fprintf('  %7.3f %10.5f  %10.2e \n',[xv y feil]');
