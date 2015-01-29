% Program tri321
% Beregner eksemplet i avsnitt 3.2.1
% , tilfelle I, med bruk av tri.
clear
beta = 2;
h = 0.1;
n = round(1.0/h); % Antall ligninger
a = ones(n,1); % Underdiagonal
b = -(2 + (beta*h)^2)*a; % Hovediagonal
c = a; % Overdiagonal
c(2) = 2.0;
d = zeros(n,1); % Høyre side
d(n) = -1;
x = tri(a,b,c,d);
x