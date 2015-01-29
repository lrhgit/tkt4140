% Program ex133
% 
% Bruker eksempel 13.3 i Ascher
clear; close;
m  = 16; l = m/2;
x = 0: pi/l : (m-1)*pi/l; % absisser på [0,2pi]
% Bruker hattefunksjonen y(x) = x for [0,pi]
% og y = 2*pi - x for [pi, 2*pi]
for j = 1:m
    fak = x(j);
    y(j) = fak;
    if fak > pi
        y(j) = 2*pi - fak;
    end
end
% Beregner koeffisientene a0 ,a og b;
[a0,a,b] = fcndft(y);
% Lager oss et fint nett på [0,2*pi]
% for interpolering.
xi = 0 : 0.02*pi : 2*pi;
% Beregner de interpolerte y-verdiene med bruk av fcndft2
yi = fcndft2(xi,a0,a,b);
% Eksakte verdier
for j = 1:length(xi)
    fak = xi(j);
    yex(j) = fak;
    if fak > pi
        yex(j) = 2*pi - fak;
    end
end
plot(xi,yex,'b--',xi,yi)
% axis([-1 2 -2 6])