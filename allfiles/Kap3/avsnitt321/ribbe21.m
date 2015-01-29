% Program ribbe21
% Beregner eksemplet i avsnitt 3.2.1,
% tilfelle II med sentraldifferanser
% med bruk av tdma.
clear
beta = 2;
h = 0.1;
n = round(1.0/h); % Antall ligninger
a = ones(n,1); % Underdiagonal
b = -(2 + (beta*h)^2)*a; % Hovediagonal
c = a; % Overdiagonal
c(1) = 2.0;
d = zeros(n,1); % Høyre side
d(n) = -1;
thet = tdma(a,b,c,d);
thet = [thet; 1];
% analytisk løsning
xa = (0 : h: 1.0)';
theta = cosh(beta*xa)/cosh(beta);
fprintf('  x      theta     analyt \n');
fprintf('%4.2f  %9.5f  %9.5f \n',[xa thet theta]');
