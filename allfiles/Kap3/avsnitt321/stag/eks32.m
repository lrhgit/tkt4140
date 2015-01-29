%============ Program Eks32 ====================
% Løser varmeledningsproblemet i avsnitt 3.2
% der diff.ligningen er gitt i ligning 3.2.4
% Bruker tdma som løser. Bruker forskjøvet nett
%===============================================
clear
h = 0.05; % skrittlengde
beta = 2.0; % Biottallet er rota av beta
n = round(1.0/h); % Antall ligninger
fact = -(2.0 + (beta*h)^2);
a = ones(n,1) ; % underdiagonal
b = fact*ones(n,1); % hoveddiagonal
c = a ; % overdiagonal 
b(1) = -(1.0 + (beta*h)^2); % for tilfelle I
b(n) = -(3.0 + (beta*h)^2); % 
d = zeros(n,1) ; % høyre side av ligningsystemet
d(n) = -2.0;
%---------------------------------------------
% Tilfelle I : lign 3.2.11 med forskjøvet nett
%---------------------------------------------
theta = tdma(a,b,c,d); % Løs ligningsystemet
%theta(n+1)= 1.0;
x = [h/2:h:1.0-h/2]';
antheta = cosh(beta*x)/cosh(beta);
rerr = (antheta - theta)./antheta; % Relativ feil
fprintf(' TILFELLE I \n\n')
fprintf('    %6.3f  %12.5f %12.3e \n',[x theta rerr]');
% %---------------------------
% % Tilfelle II : lign 3.2.14
% %---------------------------
% b(1) = 2.0;
% c(1) = -( 2.0 - (beta*h)^2);
% theta = tdma(a,b,c,d); % Løs ligningsystemet
% theta(n+1)= 1.0;
% x = [0.0:h:1.0]';
% fprintf('\n TILFELLE II \n\n')
% fprintf('    %6.2f  %12.5f \n',[x theta]');
