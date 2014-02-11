% Program ex133dct
% 
% Bruker eksempel 13.3 i Ascher, men nå med bruk av
% cosinus-transformasjon fcndct istedenfor fcndft.
clear; close;
% Bruker hattefunksjonen y(x) = x for [0,pi]
m  = 16;
iv = 0 : m -1;
x = pi*(iv + 1/2)/m;; % x-verdier på [0,pi]
y = x;
% Beregner koeffisientene a0 og a;
[a0,a] = fcndct(y);
% Lager oss et fint nett på [0,pi] for interpolering.
xi = 0 : 0.02*pi : pi;
% Beregner de interpolerte y-verdiene med bruk av fcndft2
yi = fcndct2(xi,a0,a);
% Eksakte verdier
yex = xi;
plot(xi,yex,'b--',xi,yi)
% axis([-1 2 -2 6])