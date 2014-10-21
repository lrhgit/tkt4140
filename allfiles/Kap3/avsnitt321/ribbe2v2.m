%============ Program ribbe2v2 =================
% Løser varmeledningsproblemet i avsnitt 3.2
% der diff.ligningen er gitt i ligning 3.2.4
% Som ribbe2, men bruker Matlabs innebygde
% løser for glisne matriser
%===============================================
clear
h = 0.1; % skrittlengde
beta = 2.0; % Biottallet er beta^2
n = round(1.0/h); % Antall ligninger
fact = -(2.0 + (beta*h)^2);
a = ones(n,1) ; % underdiagonal
b = fact*ones(n,1); % hoveddiagonal
c = a ; % overdiagonal 
c(2) = 2.0 ; % for versjon 1
d = zeros(n,1) ; % høyre side av ligningsystemet
d(n) = -1.0;
%----------------------------------------------
% Genererer glissen koeffisientmatrise
% for tilfelle 2 , versjon 1 : lign 2.11 + 2.13
%----------------------------------------------
A = spdiags([a b c], [-1 0 1],n,n);
theta = A\d; % Løs ligningsystemet
theta(n+1)= 1.0;
x = [0.0:h:1.0]';
fprintf(' Tilfelle 2 \n\n')
fprintf('    %6.2f  %10.5f \n',[x theta]');
%-----------------------------------
% Tilfelle 2 , versjon 2 : lign 2.17
%-----------------------------------
b(1) = 2.0;
c(2) = -( 2.0 - (beta*h)^2);
A = spdiags([a b c], [-1 0 1],n,n);
theta = A\d; % Løs ligningsystemet
theta(n+1)= 1.0;
x = [0.0:h:1.0]';
fprintf('\n Tilfelle 2 \n\n')
fprintf('    %6.2f  %10.5f \n',[x theta]');
