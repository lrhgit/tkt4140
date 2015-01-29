function x = tdma(a,b,c,d)
% Løsning av en tridiagonal matrise 
% ved bruk av Thomas-algoritmen
% Antall ligninger finnes av lengden
% av hovediagonalen b
n = length(b);
x = zeros(size(b));
% === Eliminering ===
for k = 2:n
   q = a(k)/b(k-1);
   b(k) = b(k) - c(k-1)*q;
   d(k) = d(k) - d(k-1)*q;
end
% === Innsetting ===
q = d(n)/b(n);
x(n) = q;
for k = n-1 :-1 :1
   q = (d(k) - c(k)*q)/b(k);
   x(k) = q;
end
