%======= Program eks34- deltaligninger ======
% Løser eksemplet i avsnitt 3.4
% der diff.ligningen er gitt i ligning 3.4.1
% bruker Matlab-løser
%============================================
clear
h = 0.1; % skrittlengde
n = round(1.0/h)-1; % Antall ligninger
fac = 3.0*h*h;
nitr = 6; % Antall iterasjoner
%--- a og c er konstante vektorer
a = ones(n,1) ; % underdiagonal
c = a ; % overdiagonal 
%-- b og d forandrer verdi
y = zeros(n,1); %startverdier
b = y; d = y; dy = y; 
fprintf('        Itr.      max. resid  \n');
for k = 1:nitr  
   b = -(2.0 + fac*y); % hoveddiagonal
   d = (fac*0.5)*y.^2 + 2*y; % høyre side 
   for l = 2:n-1
      d(l) = d(l) -(y(l-1) + y(l+1));
   end
   d(n) = d(n)- (1.0 + y(n-1));
   d(1) = d(1) - (4.0 + y(2));
   %--------------------------------------
   % Genererer glissen koeffisientmatrise
   %--------------------------------------
   A = spdiags([a b c], [-1 0 1],n,n);
   dy = A\d; % Løs ligningsystemet
   %--- Beregn relativ avvik
   y = y + dy;
   dymax = max(abs(dy./y));
   fprintf(' %10d     %12.3e \n',k,dymax);
 end
 %---- Utskrift av y og dy ----
 xv = [h:h:1.0-h]';
 ex = 4.0./(1 + xv).^2;
 feil = abs((y - ex)./ex);
fprintf('\n      x         y        feil  \n\n')
fprintf('  %7.3f %10.5f  %10.2e \n',[xv y feil]');
